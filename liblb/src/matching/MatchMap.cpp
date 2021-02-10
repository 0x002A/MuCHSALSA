#include "matching/MatchMap.h"

#include <algorithm>
#include <any>
#include <iterator>

#include "Util.h"
#include "graph/Graph.h"
#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

// Constants
constexpr std::size_t TH_OVERLAP = 100;

namespace lazybastard::matching {

void MatchMap::addVertexMatch(const std::string &nanoporeID, const std::string &illuminaID,
                              std::shared_ptr<VertexMatch> &&spMatch) {
  lazybastard::util::check_pointers(spMatch.get());
  std::scoped_lock<std::mutex> lck(m_vertexMutex);

  /*
   * SIDENOTE: insert will only insert if NOT present, but it will always
   * provide us with the right iterator this saves us a lot of extra checks here
   */

  // Map representing the scaffolds
  auto nanoporeIDs =
      m_scaffolds.insert(std::make_pair(illuminaID, std::map<std::string, std::shared_ptr<VertexMatch>>())).first;
  nanoporeIDs->second.insert(std::make_pair(nanoporeID, spMatch));

  // Actual map containing the matches
  auto illuminaIDs =
      m_vertexMatches.insert(std::make_pair(nanoporeID, um<std::string, std::shared_ptr<VertexMatch>>())).first;
  illuminaIDs->second.insert(std::make_pair(illuminaID, std::move(spMatch)));
}

void MatchMap::addEdgeMatch(std::string &&edgeID, const std::string &illuminaID, std::shared_ptr<EdgeMatch> &&spMatch) {
  lazybastard::util::check_pointers(spMatch.get());
  std::scoped_lock<std::mutex> lck(m_edgeMutex);

  /*
   * SIDENOTE: insert will only insert if NOT present, but it will always
   * provide us with the right iterator this saves us a lot of extra checks here
   */

  // Actual map containing the matches
  auto illuminaIDs =
      m_edgeMatches.insert(std::make_pair(std::move(edgeID), um<std::string, std::shared_ptr<EdgeMatch>>())).first;
  illuminaIDs->second.insert(std::make_pair(illuminaID, std::move(spMatch)));
}

void MatchMap::calculateEdges() {
  threading::WaitGroup wg;
  auto const jobFn = [this](threading::Job const *const pJob) { processScaffold(pJob); };

  for (auto const &[illuminaID, scaffold] : m_scaffolds) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, illuminaID, scaffold);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void MatchMap::processScaffold(gsl::not_null<threading::Job const *> const pJob) {
  auto const scaffold = std::any_cast<std::map<std::string, std::shared_ptr<VertexMatch>>>(pJob->getParam(2));

  for (auto outerIter = std::begin(scaffold); outerIter != std::end(scaffold); ++outerIter) {
    auto const *const outerMatch = outerIter->second.get();
    for (auto innerIter = std::next(outerIter, 1); innerIter != std::end(scaffold); ++innerIter) {
      auto const *const innerMatch = innerIter->second.get();
      auto const overlap = std::make_pair(std::max(outerMatch->illuminaRange.first, innerMatch->illuminaRange.first),
                                          std::min(outerMatch->illuminaRange.second, innerMatch->illuminaRange.second));

      if (overlap.first <= overlap.second && overlap.second - overlap.first > static_cast<int>(TH_OVERLAP)) {
        auto const direction = outerMatch->direction == innerMatch->direction;
        auto const isPrimary = outerMatch->isPrimary == innerMatch->isPrimary;
        auto const outerLength =
            static_cast<float>(outerMatch->illuminaRange.second - outerMatch->illuminaRange.first + 1);
        auto const innerLength =
            static_cast<float>(innerMatch->illuminaRange.second - innerMatch->illuminaRange.first + 1);
        auto const commonLength = static_cast<float>(overlap.second - overlap.first + 1);
        auto const outerScore = static_cast<float>(outerMatch->score) * (commonLength / outerLength);
        auto const innerScore = static_cast<float>(innerMatch->score) * (commonLength / innerLength);
        auto const sumScore = outerScore + innerScore;

        auto vertexIDs = std::make_pair(innerIter->first, outerIter->first);
        auto edgeID = m_pGraph->addEdge(vertexIDs);

        addEdgeMatch(std::move(edgeID), std::any_cast<std::string>(pJob->getParam(1)),
                     lazybastard::util::make_shared_aggregate<EdgeMatch>(overlap, direction, sumScore, isPrimary));
      }
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup *const>(pJob->getParam(0))->done();
}

} // namespace lazybastard::matching

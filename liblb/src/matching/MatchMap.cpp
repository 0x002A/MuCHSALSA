#include "matching/MatchMap.h"

#include <algorithm>
#include <iterator>

#include "Util.h"
#include "graph/Edge.h"
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
  const auto jobFn = [this](const threading::Job *pJob) { processScaffold(pJob); };

  for (const auto &[illuminaID, scaffold] : m_scaffolds) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, illuminaID, scaffold);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void MatchMap::processScaffold(gsl::not_null<const threading::Job *> pJob) {
  const auto scaffold = std::any_cast<std::map<std::string, std::shared_ptr<VertexMatch>>>(pJob->getParam(2));

  for (auto outerIter = scaffold.begin(); outerIter != scaffold.end(); ++outerIter) {
    auto *const outerMatch = outerIter->second.get();
    for (auto innerIter = std::next(outerIter, 1); innerIter != scaffold.end(); ++innerIter) {
      auto *const innerMatch = innerIter->second.get();
      const auto overlap = std::make_pair(std::max(outerMatch->illuminaRange.first, innerMatch->illuminaRange.first),
                                          std::min(outerMatch->illuminaRange.second, innerMatch->illuminaRange.second));

      if (overlap.first <= overlap.second && overlap.second - overlap.first > TH_OVERLAP) {
        const auto direction = outerMatch->direction == innerMatch->direction;
        const auto isPrimary = outerMatch->isPrimary == innerMatch->isPrimary;
        const auto outerLength = outerMatch->illuminaRange.second - outerMatch->illuminaRange.first + 1;
        const auto innerLength = innerMatch->illuminaRange.second - innerMatch->illuminaRange.first + 1;
        const auto commonLength = overlap.second - overlap.first + 1;
        const auto outerScore = outerMatch->score * commonLength / outerLength;
        const auto innerScore = innerMatch->score * commonLength / innerLength;
        const auto sumScore = outerScore + innerScore;

        auto vertexIDs = std::make_pair(innerIter->first, outerIter->first);
        auto edgeID = m_pGraph->addEdge(vertexIDs);

        addEdgeMatch(std::move(edgeID), std::any_cast<std::string>(pJob->getParam(1)),
                     lazybastard::util::make_shared_aggregate<EdgeMatch>(overlap, direction, sumScore, isPrimary));
      }
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup *>(pJob->getParam(0))->done();
}

} // namespace lazybastard::matching

#include "matching/MatchMap.h"

#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"
#include "threading/Job.h"
#include "graph/Graph.h"
#include "graph/Edge.h"

// Constants
constexpr std::size_t TH_OVERLAP = 100;

namespace lazybastard {
namespace matching {

void
MatchMap::addVertexMatch(
            const std::string& nanoporeID,
            const std::string& illuminaID,
            std::shared_ptr<VertexMatch>&& match)
{
  std::scoped_lock<std::mutex> lck(m_vertexMutex);

  /*
   * SIDENOTE: insert will only insert if NOT present, but it will always provide us with the right iterator
   *           this saves us a lot of extra checks here
   */

   // Map representing the scaffolds
   auto nanoporeIDs = m_scaffolds.insert(std::make_pair(illuminaID, std::map<std::string, std::shared_ptr<VertexMatch>>())).first;
   nanoporeIDs->second.insert(std::make_pair(nanoporeID, match));

   // Actual map containing the matches
   auto illuminaIDs = m_vertexMatches.insert(std::make_pair(nanoporeID, um<std::string, std::shared_ptr<VertexMatch>>())).first;
   illuminaIDs->second.insert(std::make_pair(illuminaID, std::move(match)));
}

void
MatchMap::addEdgeMatch(const std::string& edgeID, const std::string& illuminaID, EdgeMatch&& match)
{
  std::scoped_lock<std::mutex> lck(m_edgeMutex);

  /*
   * SIDENOTE: insert will only insert if NOT present, but it will always provide us with the right iterator
   *           this saves us a lot of extra checks here
   */

   // Actual map containing the matches
   auto illuminaIDs = m_edgeMatches.insert(std::make_pair(edgeID, um<std::string, EdgeMatch>())).first;
   illuminaIDs->second.insert(std::make_pair(illuminaID, std::move(match)));
}

void
MatchMap::calculateEgdes(gsl::not_null<threading::ThreadPool*> pThreadPool,
                         gsl::not_null<graph::Graph*> pGraph)
{
  threading::WaitGroup wg;
  const auto jobFn = [this](const threading::Job* pJob) { processScaffold(pJob); };

  for (const auto& [illuminaID, scaffold] : m_scaffolds)
  {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, illuminaID, scaffold);
    pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void
MatchMap::processScaffold(gsl::not_null<const threading::Job*> pJob)
{
  const auto scaffold = std::any_cast<std::map<std::string, std::shared_ptr<VertexMatch>>>(pJob->getParam(2));

  for (auto outerIter = scaffold.begin(); outerIter != scaffold.end(); ++outerIter)
  {
    const auto outerMatch = outerIter->second.get();
    for (auto innerIter = std::next(outerIter, 1); innerIter != scaffold.end(); ++innerIter)
    {
      const auto innerMatch = innerIter->second.get();
      const auto overlap = std::make_pair(std::max(outerMatch->illuminaRange.first,
                                                   innerMatch->illuminaRange.first),
                                          std::min(outerMatch->illuminaRange.second,
                                                   innerMatch->illuminaRange.second));

      if (overlap.first <= overlap.second && overlap.second - overlap.first > TH_OVERLAP) {
        const auto direction = outerMatch->direction == innerMatch->direction;
        const auto thresholdsPassed = outerMatch->thresholdsPassed == innerMatch->thresholdsPassed;
        const auto outerLength = outerMatch->illuminaRange.second - outerMatch->illuminaRange.first + 1;
        const auto innerLength = innerMatch->illuminaRange.second - innerMatch->illuminaRange.first + 1;
        const auto commonLength = overlap.second - overlap.first + 1;
        const auto outerScore = outerMatch->score * commonLength/outerLength;
        const auto innerScore = innerMatch->score * commonLength/innerLength;
        const auto minScore = std::min(outerScore, innerScore);
        const auto sumScore = outerScore + innerScore;

        m_pGraph->addEdge(std::make_pair(innerIter->first, outerIter->first));

        addEdgeMatch(lazybastard::graph::Edge::getEdgeID(innerIter->first, outerIter->first),
                     std::any_cast<std::string>(pJob->getParam(1)),
                     {overlap, direction, sumScore, thresholdsPassed});
      }
    }
  }

  std::any_cast<lazybastard::threading::WaitGroup*>(pJob->getParam(0))->done();
}

}
}

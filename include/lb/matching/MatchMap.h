#pragma once

#include <utility>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <map>
#include <mutex>
#include <gsl/pointers>

#include <atomic>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace threading {
  class ThreadPool;
  class Job;
}
namespace graph {
  class Graph;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

namespace matching {

/**
 * Struct representing a vertex match.
 */
struct VertexMatch {
  VertexMatch(std::pair<int, int> nanoporeRange, std::pair<int, int> illuminaRange, float rRatio, bool direction,
              std::size_t score, bool thresholdsPassed)
  : nanoporeRange(nanoporeRange)
  , illuminaRange(illuminaRange)
  , rRatio(rRatio)
  , direction(direction)
  , score(score)
  , thresholdsPassed(thresholdsPassed) {};
  const std::pair<int, int> nanoporeRange; /*!< Nanopore range */
  const std::pair<int, int> illuminaRange; /*!< Illumina range */
  const float rRatio; /*!< Read ratio */
  const bool direction; /*!< Read direction */
  const std::size_t score; /*!< Score (number of matches) */
  const bool thresholdsPassed; /*!< Did the match pass the thresholds concering length and match count */
};

/**
 * Struct representing an edge match.
 */
struct EdgeMatch {
  const std::pair<int, int> overlap; /*!< Overlap */
  const bool direction; /*!< Edge direction */
  const std::size_t score; /*!< Score */
  const bool thresholdsPassed; /*!< Did the match pass the thresholds */
};

/**
 * Class representing a map containing illumina ids matched with nanopore ids.
 *
 * Instances of this class are designed to be thread-safe.
 */
class MatchMap {
public:
  MatchMap(lazybastard::graph::Graph* pGraph)
    :m_pGraph(pGraph) {};
  /**
   * Adds a vertex match to the map.
   *
   * @param nanoporeID the nanopore id
   * @param illuminaID the illumina id
   * @param match the node match to be added to the map
   */
  void addVertexMatch(
        const std::string& nanoporeID,
        const std::string& illuminaID,
        std::shared_ptr<VertexMatch>&& match);

  /**
   * Adds an edge match to the map.
   *
   * @param edgeID the identifier of the edge
   * @param illuminaID the illumina id
   * @param match the edge match to be added to the map
   */
  void addEdgeMatch(const std::string& edgeID, const std::string& illuminaID, EdgeMatch&& match);

  /**
   * Calculates or updates Edges according to the scaffolds.
   *
   * @param pThreadPool pointer to the thread pool to be used for parallelization
   * @param pGraph pointer to the graph receiving the edges
   */
  void calculateEgdes(
        gsl::not_null<lazybastard::threading::ThreadPool*> pThreadPool,
        gsl::not_null<lazybastard::graph::Graph*> pGraph);

  /**
   * Processes a scaffold.
   *
   * @param pJob pointer to the job containing the parameters
   */
  void processScaffold(gsl::not_null<const threading::Job*> pJob);
private:
  template<typename T1, typename T2> using um = std::unordered_map<T1, T2>;

  um<std::string, um<std::string, std::shared_ptr<VertexMatch>>> m_vertexMatches; /*!< Map containing the node matches */
  um<std::string, um<std::string, EdgeMatch>> m_edgeMatches; /*!< Map containing the edge matches */
  um<std::string, std::map<std::string, std::shared_ptr<VertexMatch>>> m_scaffolds; /*!< Map containing the scaffolds */
  std::mutex m_vertexMutex; /*!< Mutex for securing the parallel use of the edge match map */
  std::mutex m_edgeMutex; /*!< Mutex for securing the parallel use of the edge match map */

  lazybastard::graph::Graph* m_pGraph;
};

}
}

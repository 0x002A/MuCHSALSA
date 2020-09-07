#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace threading {
class ThreadPool;
class Job;
} // namespace threading
namespace graph {
class Graph;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

namespace matching {

/**
 * Struct representing a match attached to a Vertex.
 */
struct VertexMatch {
  const std::pair<int, int> nanoporeRange; /*!< Nanopore range */
  const std::pair<int, int> illuminaRange; /*!< Illumina range */
  const float rRatio;                      /*!< Read ratio */
  const bool direction;                    /*!< Read direction */
  const std::size_t score;                 /*!< Score (number of matches) */
  const bool thresholdsPassed;             /*!< Did the match pass the thresholds concering
                                              length and match count */
};

/**
 * Struct representing a match attached to an Edge.
 */
struct EdgeMatch {
  const std::pair<int, int> overlap; /*!< Overlap */
  const bool direction;              /*!< Edge direction */
  const std::size_t score;           /*!< Score */
  const bool thresholdsPassed;       /*!< Did the match pass the thresholds */
};

/**
 * Class representing a map containing illumina IDs matched with nanopore IDs.
 *
 * Instances of this class are designed to be thread-safe.
 */
class MatchMap {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @param pThreadPool pointer to the ThreadPool to be used for parallelization
   * @param pGraph pointer to the Graph receiving the Vertex instances
   */
  MatchMap(gsl::not_null<threading::ThreadPool *> pThreadPool, gsl::not_null<graph::Graph *> pGraph)
      : m_pThreadPool(pThreadPool), m_pGraph(pGraph){};
  /**
   * Adds a Vertex match to the map.
   *
   * @param nanoporeID the nanopore ID
   * @param illuminaID the illumina ID
   * @param spMatch shared pointer to the Vertex match to be added to the map
   */
  void addVertexMatch(const std::string &nanoporeID, const std::string &illuminaID,
                      std::shared_ptr<VertexMatch> &&spMatch);

  /**
   * Adds an Edge match to the map.
   *
   * @param edgeID the identifier of the Edge
   * @param illuminaID the illumina ID
   * @param spMatch shared pointer to the Edge match to be added to the map
   */
  void addEdgeMatch(const std::string &edgeID, const std::string &illuminaID, std::shared_ptr<EdgeMatch> &&spMatch);

  /**
   * Creates or updates Edge instances according to the scaffolds.
   */
  void calculateEgdes();

  /**
   * Processes a scaffold.
   *
   * @param pJob pointer to the job containing the parameters
   */
  void processScaffold(gsl::not_null<const threading::Job *> pJob);

private:
  template <typename T1, typename T2> using um = std::unordered_map<T1, T2>;

  um<std::string, um<std::string, std::shared_ptr<VertexMatch>>>
      m_vertexMatches;                                                        /*!< Map containing the Vertex matches */
  um<std::string, um<std::string, std::shared_ptr<EdgeMatch>>> m_edgeMatches; /*!< Map containing the Edge matches */
  um<std::string, std::map<std::string, std::shared_ptr<VertexMatch>>> m_scaffolds; /*!< Map containing the scaffolds */
  std::mutex m_vertexMutex;             /*!< Mutex for securing the parallel use of the Edge MatchMap */
  std::mutex m_edgeMutex;               /*!< Mutex for securing the parallel use of the Edge MatchMap */
  threading::ThreadPool *m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  graph::Graph *m_pGraph;               /*!< Pointer to the Graph receiving the vertices */
};

} // namespace matching
} // namespace lazybastard

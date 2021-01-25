#pragma once

#include <gsl/pointers>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace threading {
class ThreadPool;
class Job;
} // namespace threading
namespace graph {
class Graph;
} // namespace graph
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

namespace matching {

/**
 * Struct representing a match attached to a Vertex.
 */
struct VertexMatch {
  std::pair<int, int> const nanoporeRange; /*!< Nanopore range */
  std::pair<int, int> const illuminaRange; /*!< Illumina range */
  float const rRatio;                      /*!< Read ratio */
  bool const direction;                    /*!< Read direction */
  std::size_t const score;                 /*!< Score (number of matches) */
  bool const isPrimary;                    /*!< Did the match pass the thresholds concerning
                                                     length and match count */
};

/**
 * Struct representing a match attached to an Edge.
 */
struct EdgeMatch {
  std::pair<int, int> const overlap; /*!< Overlap */
  bool const direction;              /*!< Edge direction */
  std::size_t const score;           /*!< Score */
  bool const isPrimary;              /*!< Did the match pass the thresholds */
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
  MatchMap(gsl::not_null<threading::ThreadPool *> const pThreadPool, gsl::not_null<graph::Graph *> const pGraph)
      : m_pThreadPool(pThreadPool), m_pGraph(pGraph){};
  /**
   * Adds a Vertex match to the map.
   *
   * @param nanoporeID the nanopore ID
   * @param illuminaID the illumina ID
   * @param spMatch shared pointer to the Vertex match to be added to the map
   */
  void addVertexMatch(std::string const &nanoporeID, std::string const &illuminaID,
                      std::shared_ptr<VertexMatch> &&spMatch);

  /**
   * Getter for a specific match of the Vertex having the supplied ID.
   *
   * @param vertexID the identifier of the Vertex
   * @param illuminaID the illumina ID to return the match for
   * @return A pointer to the VertexMatch if found nullptr otherwise
   */
  [[nodiscard]] VertexMatch const *getVertexMatch(std::string const &vertexID, std::string const &illuminaID) const {
    auto const &vertexIter = m_vertexMatches.find(vertexID);
    if (vertexIter != m_vertexMatches.end()) {
      auto const &illuminaIter = vertexIter->second.find(illuminaID);
      if (illuminaIter != vertexIter->second.end()) {
        return illuminaIter->second.get();
      }
    }

    return nullptr;
  };

  /**
   * Adds an Edge match to the map.
   *
   * @param edgeID the identifier of the Edge
   * @param illuminaID the illumina ID
   * @param spMatch shared pointer to the Edge match to be added to the map
   */
  void addEdgeMatch(std::string &&edgeID, std::string const &illuminaID, std::shared_ptr<EdgeMatch> &&spMatch);

  /**
   * Getter for Edge matches.
   *
   * @return The map containing the Edge matches
   */
  [[nodiscard]] auto const &getEdgeMatches() const { return m_edgeMatches; };

  /**
   * Creates or updates Edge instances according to the scaffolds.
   */
  void calculateEdges();

  /**
   * Processes a scaffold.
   *
   * @param pJob pointer to the job containing the parameters
   */
  void processScaffold(gsl::not_null<threading::Job const *> pJob);

private:
  template <typename T1, typename T2> using um = std::unordered_map<T1, T2>;

  um<std::string, um<std::string, std::shared_ptr<VertexMatch>>>
      m_vertexMatches;                                                        /*!< Map containing the Vertex matches */
  um<std::string, um<std::string, std::shared_ptr<EdgeMatch>>> m_edgeMatches; /*!< Map containing the Edge matches */
  um<std::string, std::map<std::string, std::shared_ptr<VertexMatch>>> m_scaffolds; /*!< Map containing the scaffolds */
  std::mutex m_vertexMutex;                   /*!< Mutex for securing the parallel use of the Edge MatchMap */
  std::mutex m_edgeMutex;                     /*!< Mutex for securing the parallel use of the Edge MatchMap */
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  graph::Graph *const m_pGraph;               /*!< Pointer to the Graph receiving the vertices */
};

} // namespace matching
} // namespace lazybastard

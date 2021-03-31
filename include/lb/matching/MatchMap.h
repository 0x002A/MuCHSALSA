#pragma once

#include <cstddef>
#include <gsl/pointers>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>

#include "Lb.fwd.h"
#include "types/Toggle.h"

namespace lazybastard {

//// UTILITY FUNCTIONS ////
struct MatchingUtil {
  static bool scaffoldLineIdxCmp(graph::Vertex const *pV1, graph::Vertex const *pV2);
};
//// /UTILITY FUNCTIONS ////

namespace matching {

/**
 * Struct representing a match attached to a Vertex.
 */
struct VertexMatch {
  std::pair<int, int> const nanoporeRange; /*!< Nanopore range */
  std::pair<int, int> const illuminaRange; /*!< Illumina range */
  double const rRatio;                     /*!< Read ratio */
  Toggle const direction;                  /*!< Read direction */
  std::size_t const score;                 /*!< Score (number of matches) */
  Toggle const isPrimary;                  /*!< Did the match pass the thresholds concerning
                                              length and match count */
};

/**
 * Struct representing a match attached to an Edge.
 */
struct EdgeMatch {
  std::pair<int, int> const overlap; /*!< Overlap */
  Toggle const direction;            /*!< Edge direction */
  double const score;                /*!< Score */
  Toggle const isPrimary;            /*!< Did the match pass the thresholds */
};

/**
 * Class representing a map containing Illumina IDs matched with Nanopore IDs.
 *
 * Instances of this class are designed to be thread-safe.
 */
class MatchMap {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param pThreadPool a pointer to the ThreadPool to be used for parallelization
   * @param pGraph a pointer to the Graph receiving the Vertex instances
   */
  MatchMap(gsl::not_null<threading::ThreadPool *> const pThreadPool, gsl::not_null<graph::Graph *> const pGraph)
      : m_pThreadPool(pThreadPool), m_pGraph(pGraph){};
  /**
   * Adds a VertexMatch to the map.
   *
   * @param nanoporeID a constant reference to the nanopore ID
   * @param illuminaID a constant reference to the illumina ID
   * @param spMatch an rvalue reference to the std::shared_ptr to the VertexMatch to be added to the map (by moving)
   */
  void addVertexMatch(std::string const &nanoporeID, std::string const &illuminaID,
                      std::shared_ptr<VertexMatch> &&spMatch);

  /**
   * Getter returning a specific match of the Vertex having the supplied ID.
   * This functions returns nullptr if the requested VertexMatch wasn't found.
   *
   * @param vertexID a constant reference to the identifier of the Vertex
   * @param illuminaID a constant reference to the illumina ID associated with the VertexMatch
   * @return A pointer to the VertexMatch (constant) if found nullptr otherwise
   */
  [[nodiscard]] VertexMatch const *getVertexMatch(std::string const &vertexID, std::string const &illuminaID) const {
    auto const vertexIter = m_vertexMatches.find(vertexID);
    if (vertexIter != std::end(m_vertexMatches)) {
      auto const illuminaIter = vertexIter->second.find(illuminaID);
      if (illuminaIter != std::end(vertexIter->second)) {
        return illuminaIter->second.get();
      }
    }

    return nullptr;
  };

  /**
   * Getter returning all EdgeMatch instances for a specific Edge stored within this map.
   *
   * @param vertexID a constant reference to the identifier of the Vertex
   * @return A pointer to the std::unordered_map containing all VertexMatch instances associated with the requested
   *         Vertex
   */
  [[nodiscard]] std::unordered_map<std::string, std::shared_ptr<VertexMatch>> const *
  getVertexMatches(std::string const &vertexID) const {
    auto iter = m_vertexMatches.find(vertexID);
    if (iter != std::end(m_vertexMatches)) {
      return &(iter->second);
    }

    return nullptr;
  };

  /**
   * Getter returning all VertexMatch instances stored within this map.
   *
   * @return The std::unordered_map containing all VertexMatch instances
   */
  [[nodiscard]] auto const &getVertexMatches() const { return m_vertexMatches; };

  /**
   * Adds an Edge match to the map.
   *
   * @param edgeID an rvalue reference to the identifier of the Edge
   * @param illuminaID a constant reference the illumina ID
   * @param spMatch an rvalue reference to the std::shared_pointer to the EdgeMatch to be added to the map (by moving)
   */
  void addEdgeMatch(std::string &&edgeID, std::string const &illuminaID, std::shared_ptr<EdgeMatch> &&spMatch);

  /**
   * Getter returning a specific match of the Edge having the supplied ID.
   * This functions returns nullptr if the requested EdgeMatch wasn't found.
   *
   * @param edgeID a constant reference to the identifier of the Edge
   * @param illuminaID a constant reference to the illumina ID associated with the EdgeMatch
   * @return A pointer to the EdgeMatch (constant) if found nullptr otherwise
   */
  [[nodiscard]] EdgeMatch const *getEdgeMatch(std::string const &edgeID, std::string const &illuminaID) const {
    auto const edgeIter = m_edgeMatches.find(edgeID);
    if (edgeIter != std::end(m_edgeMatches)) {
      auto const illuminaIter = edgeIter->second.find(illuminaID);
      if (illuminaIter != std::end(edgeIter->second)) {
        return illuminaIter->second.get();
      }
    }

    return nullptr;
  };

  /**
   * Getter returning all EdgeMatch instances for a specific Edge stored within this map.
   *
   * @param edgeID a constant reference to the identifier of the Edge
   * @return A pointer to the std::unordered_map containing all EdgeMatch instances associated with the requested Edge
   */
  [[nodiscard]] std::unordered_map<std::string, std::shared_ptr<EdgeMatch>> const *
  getEdgeMatches(std::string const &edgeID) const {
    auto const iter = m_edgeMatches.find(edgeID);
    if (iter != std::end(m_edgeMatches)) {
      return &iter->second;
    }

    return nullptr;
  };

  /**
   * Getter returning all EdgeMatch instances stored within this map.
   *
   * @return The std::unordered_map containing all EdgeMatch instances
   */
  [[nodiscard]] auto const &getEdgeMatches() const { return m_edgeMatches; };

  /**
   * Creates or updates Edge instances according to the scaffolds.
   */
  void calculateEdges();

  /**
   * Processes a scaffold.
   *
   * @param pJob a pointer to the threading::Job containing the parameters
   */
  void processScaffold(gsl::not_null<threading::Job const *> pJob);

private:
  template <typename T1, typename T2> using um = std::unordered_map<T1, T2>;

  um<std::string, um<std::string, std::shared_ptr<VertexMatch>>>
      m_vertexMatches; /*!< Map containing the VertexMatch instances */
  um<std::string, um<std::string, std::shared_ptr<EdgeMatch>>>
      m_edgeMatches; /*!< Map containing the EdgeMatch instances */
  um<std::string,
     std::map<graph::Vertex *, std::shared_ptr<VertexMatch>, decltype(&lazybastard::MatchingUtil::scaffoldLineIdxCmp)>>
      m_scaffolds;          /*!< Map containing the scaffolds */
  std::mutex m_vertexMutex; /*!< std::mutex for securing the parallel use of the map containing VertexMatches */
  std::mutex m_edgeMutex;   /*!< std::mutex for securing the parallel use of the map containing EdgeMatches */
  threading::ThreadPool *const m_pThreadPool; /*!< Pointer to the ThreadPool used for parallelization */
  graph::Graph *const m_pGraph;               /*!< Pointer to the Graph receiving the Vertex instances */
};

} // namespace matching
} // namespace lazybastard

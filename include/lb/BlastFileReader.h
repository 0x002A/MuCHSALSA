#pragma once

#include <gsl/pointers>
#include <iosfwd>

#include "Lb.fwd.h"

namespace lazybastard {

/**
 * Class representing a reader for reading files formatted according to the so called BLAST format.
 *
 * The reader creates Vertex instances and adds them to the supplied Graph.
 */
class BlastFileReader {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param pThreadPool a pointer to the ThreadPool to be used for parallelization
   * @param inputStream the input stream of the file
   * @param pGraph a pointer to the Graph receiving the Vertex instances
   * @param pMatchMap a pointer to the MatchMap
   */
  BlastFileReader(gsl::not_null<threading::ThreadPool *> const pThreadPool, std::ifstream &inputStream,
                  gsl::not_null<graph::Graph *> const pGraph, gsl::not_null<matching::MatchMap *> const pMatchMap)
      : m_pThreadPool(pThreadPool), m_inputStream(inputStream), m_pGraph(pGraph), m_pMatchMap(pMatchMap){};

  /**
   * Destructor.
   */
  ~BlastFileReader() = default;

  /**
   * Copying is disallowed.
   */
  BlastFileReader(BlastFileReader const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  BlastFileReader &operator=(BlastFileReader const &) = delete;

  /**
   * Moving is disallowed.
   */
  BlastFileReader(BlastFileReader &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  BlastFileReader &operator=(BlastFileReader &&) = delete;

  /**
   * Reads the file.
   */
  void read();

  /**
   * Parses a line of the file.
   *
   * @param pJob a pointer to the Job containing the parameters
   */
  void parseLine(gsl::not_null<threading::Job const *> pJob);

private:
  std::ifstream &m_inputStream;               /*!< Input stream of the file */
  threading::ThreadPool *const m_pThreadPool; /*!< A pointer to the ThreadPool used for parallelization */
  graph::Graph *const m_pGraph;               /*!< A pointer to the Graph receiving the Vertex instances */
  matching::MatchMap *const m_pMatchMap;      /*!< A pointer to the MatchMap */
};

} // namespace lazybastard

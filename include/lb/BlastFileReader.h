#pragma once

#include <iosfwd>
#include <gsl/pointers>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace threading {
  class ThreadPool;
  class Job;
}
namespace graph {
  class Graph;
}
namespace matching {
  class MatchMap;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * Class representing a reader for reading files formatted according to the so called BLAST format.
 *
 * The reader creates vertex objects and adds them to the supplied graph.
 */
class BlastFileReader {
public:
  /**
   * Class constructor which creates a new instance.
   *
   * @param pThreadPool pointer to the thread pool to be used for parallelization
   * @param inputStream input stream of the file
   * @param pGraph pointer to the graph receiving the vertices
   * @param pMatchMap pointer to the match map
   */
  BlastFileReader(gsl::not_null<threading::ThreadPool*> pThreadPool, std::ifstream& inputStream,
                  gsl::not_null<graph::Graph*> pGraph, gsl::not_null<matching::MatchMap*> pMatchMap)
    : m_pThreadPool(pThreadPool)
    , m_inputStream(inputStream)
    , m_pGraph(pGraph)
    , m_pMatchMap(pMatchMap) {};

  /**
   * Reads the file.
   */
  void read();

  /**
   * Parses a line of the file.
   *
   * @param pJob pointer to the job containing the parameters
   */
  void parseLine(gsl::not_null<const threading::Job*> pJob);
protected:
  std::ifstream& m_inputStream; /*!< Input stream of the file */
  threading::ThreadPool* m_pThreadPool; /*!< Pointer to the thread pool used for parallelization */
  graph::Graph* m_pGraph; /*!< Pointer to the graph receiving the vertices */
  matching::MatchMap* m_pMatchMap; /*!< Pointer to the match map */
};

}

#pragma once

#include <cstddef>
#include <string>

/**
 * Class for parameters.
 *
 * This class holds the input parameters and provides application-specific functions.
 * Instances of this class are immutable by default and therefore thread-safe.
 */
class Application {
public:
  /**
   * Class constructor which creates a new instance.
   */
  Application(int argc, char *argv[]);

  /**
   * Copying is disallowed.
   */
  Application(const Application &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Application &operator=(const Application &) = delete;

  /**
   * Moving is disallowed.
   */
  Application(Application &&) = delete;

  /**
   * Move assignment is disallowed.
   */
  Application &operator=(Application &&) = delete;

  /**
   * Checks if files and directories exist and whether they are readable and
   * writeable.
   */
  bool checkIntegrity() const;

  /**
   * Getter for the supplied path to the file containing the contigs.
   *
   * @return The supplied path to the file containing the contigs
   */
  const std::string &getContigsFilePath() const noexcept { return m_contigsFilePath; };

  /**
   * Getter for the supplied path to the file containing the unitigs.
   *
   * @return The supplied path to the file containing the unitigs
   */
  const std::string &getUnitigsFilePath() const noexcept { return m_unitigsFilePath; };

  /**
   * Getter for the supplied path to the file containing the nanopore data.
   *
   * @return The supplied path to the file containing the nanopore data
   */
  const std::string &getNanoporeFilePath() const noexcept { return m_nanoporeFilePath; };

  /**
   * Getter for the supplied path to the output directory.
   *
   * @return The supplied path to the output directory
   */
  const std::string &getOutputFilePath() const noexcept { return m_outputPath; };

  /**
   * Getter for the supplied level of parallelism.
   *
   * @return The supplied level of parallelism
   */
  const std::size_t &getThreadCount() const noexcept { return m_threadCount; };

private:
  //@{
  //** Filepath */
  std::string m_contigsFilePath, m_unitigsFilePath, m_nanoporeFilePath, m_outputPath;
  //@}
  std::size_t m_threadCount; /*!< Number of threads to use */

  /**
   * Parses the parameters from the command line.
   * @param argc number of parameters.
   * @param argv array of parameters.
   */
  void readParameters(int argc, char *argv[]);
};

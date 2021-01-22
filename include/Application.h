#pragma once

#include <cstddef>
#include <gsl/span>
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
  explicit Application(gsl::span<char const *> const &args);

  /**
   * Destructor.
   */
  ~Application() = default;

  /**
   * Copying is disallowed.
   */
  Application(Application const &) = delete;

  /**
   * Copy assignment is disallowed.
   */
  Application &operator=(Application const &) = delete;

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
  [[nodiscard]] auto checkIntegrity() const -> bool;

  /**
   * Getter for the supplied path to the file containing the contigs.
   *
   * @return The supplied path to the file containing the contigs
   */
  [[nodiscard]] std::string const &getContigsFilePath() const noexcept { return m_contigsFilePath; };

  /**
   * Getter for the supplied path to the file containing the unitigs.
   *
   * @return The supplied path to the file containing the unitigs
   */
  [[nodiscard]] std::string const &getUnitigsFilePath() const noexcept { return m_unitigsFilePath; };

  /**
   * Getter for the supplied path to the file containing the nanopore data.
   *
   * @return The supplied path to the file containing the nanopore data
   */
  [[nodiscard]] std::string const &getNanoporeFilePath() const noexcept { return m_nanoporeFilePath; };

  /**
   * Getter for the supplied path to the output directory.
   *
   * @return The supplied path to the output directory
   */
  [[nodiscard]] std::string const &getOutputFilePath() const noexcept { return m_outputPath; };

  /**
   * Getter for the supplied level of parallelism.
   *
   * @return The supplied level of parallelism
   */
  [[nodiscard]] std::size_t const &getThreadCount() const noexcept { return m_threadCount; };

  /**
   * Getter for the supplied wiggle room.
   *
   * @return The supplied level of wiggle room
   */
  [[nodiscard]] std::size_t const &getWiggleRoom() const noexcept { return m_wiggleRoom; };

private:
  //@{
  //** Filepath */
  std::string m_contigsFilePath, m_unitigsFilePath, m_nanoporeFilePath, m_outputPath;
  //@}
  std::size_t m_threadCount;     /*!< Number of threads to use */
  std::size_t m_wiggleRoom{700}; /*!< Wiggle room */

  /**
   * Parses the parameters from the command line.
   * @param argc number of parameters.
   * @param argv array of parameters.
   */
  void readParameters(gsl::span<char const *> const &args);
};

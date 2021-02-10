#pragma once

#include <iosfwd>
#include <mutex>

namespace lazybastard {

/**
 * Class offering **thread-safe** writing to output streams.
 */
class OutputWriter {
public:
  /**
   * Class constructor creating a new instance.
   *
   * @param osQuery a reference to the query output stream
   * @param osPAF a reference to the PAF output stream
   * @param osTarget a reference to the target output stream
   */
  OutputWriter(std::ostream &osQuery, std::ostream &osPAF, std::ostream &osTarget);

private:
  //@{
  //** Reference to a specific output stream */
  std::ostream &m_osQuery, &m_osPAF, &m_osTarget;
  //@}
  //@{
  //** Mutex for securing the corresponding output stream */
  std::mutex m_mutexQuery, m_mutexPAF, m_mutexTarget;
  //@}
};

} // namespace lazybastard
#include "OutputWriter.h"

#include <exception>
#include <iostream>

namespace lazybastard {
OutputWriter::OutputWriter(std::ostream &osQuery, std::ostream &osPAF, std::ostream &osTarget)
    : m_osQuery(osQuery), m_osPAF(osPAF), m_osTarget(osTarget) {
  if (!m_osQuery.good() || !m_osPAF.good() || !m_osTarget.good()) {
    throw std::runtime_error("Can't write to output files! Aborting!");
  }
}
} // namespace lazybastard
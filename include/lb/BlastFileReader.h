#pragma once

#include <iosfwd>
#include <gsl/pointers>

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace threading {
class ThreadPool;
class Job;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class BlastFileReader {
public:
  BlastFileReader(std::ifstream& inputStream, gsl::not_null<threading::ThreadPool*> pThreadPool);
  void parseLine(gsl::not_null<const threading::Job*> pJob);
protected:
  std::ifstream& m_inputStream;
  threading::ThreadPool* m_pThreadPool;
};

}

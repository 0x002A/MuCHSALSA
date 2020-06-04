#include <fstream>
#include <functional>
#include <string>

#include "BlastFileReader.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"
#include "threading/Job.h"

namespace lazybastard {

BlastFileReader::BlastFileReader(std::ifstream& inputStream, gsl::not_null<threading::ThreadPool*> pThreadPool)
  : m_inputStream(inputStream)
  , m_pThreadPool(pThreadPool)
{
  threading::WaitGroup wg;
  auto jobFn = [this](const threading::Job* pJob) { parseLine(pJob); };

  std::string line;
  while (std::getline(m_inputStream, line))
  {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, line);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void
BlastFileReader::parseLine(gsl::not_null<const threading::Job*> pJob)
{
  // to be implemented

  std::any_cast<threading::WaitGroup*>(pJob->getParam(0))->done();
}

}

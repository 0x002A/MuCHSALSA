#include "BlastFileReader.h"

#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"
#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

// Constants
constexpr std::size_t MINIMUM_MATCHES = 400;
constexpr std::size_t TH_LENGTH = 500;
constexpr std::size_t TH_MATCHES = 500;

namespace lazybastard {

void BlastFileReader::read() {
  threading::WaitGroup wg;
  auto jobFn = [this](const threading::Job *pJob) { parseLine(pJob); };

  std::string line;
  while (std::getline(m_inputStream, line)) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, line);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void BlastFileReader::parseLine(gsl::not_null<const threading::Job *> pJob) {
  std::vector<std::string> tokens;

  std::istringstream iss(std::any_cast<std::string>(pJob->getParam(1)), std::ios_base::in);
  std::string token;
  while (std::getline(iss, token, '\t')) {
    tokens.push_back(token);
  }

  if (tokens.size() < 15) {
    throw std::runtime_error("Invalid BLAST file.");
  }

  const auto illuminaRange = std::make_pair(std::stoi(tokens[2]), std::stoi(tokens[3]) - 1);
  const auto matches = static_cast<std::size_t>(std::stoi(tokens[9]));

  const auto nanoporeLength = std::stoi(tokens[6]);

  auto addNode = matches >= MINIMUM_MATCHES;
  addNode &= illuminaRange.second - illuminaRange.first + 1 >= MINIMUM_MATCHES;

  if (addNode) {
    auto spVertex = std::make_shared<graph::Vertex>(tokens[5], nanoporeLength);
    m_pGraph->addVertex(std::move(spVertex));

    const auto &nanoporeID = tokens[5];
    const auto &illuminaID = tokens[0];

    const auto nanoporeRange = std::make_pair(std::stoi(tokens[7]), std::stoi(tokens[8]) - 1);
    const auto direction = tokens[4] == "+";
    const auto rRatio = static_cast<float>(illuminaRange.second - illuminaRange.first + 1) /
                        (nanoporeRange.second - nanoporeRange.first + 1);

    auto thresholdsPassed = illuminaRange.second - illuminaRange.first + 1 >= TH_LENGTH;
    thresholdsPassed &= matches >= TH_MATCHES;

    auto spVertexMatch = lazybastard::util::make_shared_aggregate<lazybastard::matching::VertexMatch>(
        nanoporeRange, illuminaRange, rRatio, direction, matches, thresholdsPassed);
    m_pMatchMap->addVertexMatch(nanoporeID, illuminaID, std::move(spVertexMatch));
  }

  std::any_cast<threading::WaitGroup *>(pJob->getParam(0))->done();
}

} // namespace lazybastard

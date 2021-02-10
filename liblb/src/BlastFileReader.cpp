#include "BlastFileReader.h"

#include <algorithm>
#include <any>
#include <cstddef>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"
#include "threading/Job.h"
#include "threading/ThreadPool.h"
#include "threading/WaitGroup.h"

// Constants
constexpr static std::size_t MINIMUM_MATCHES = 400;
constexpr static std::size_t TH_LENGTH = 500;
constexpr static std::size_t TH_MATCHES = 500;

constexpr static std::size_t POS_IID = 0;
constexpr static std::size_t POS_NID = 5;
constexpr static std::size_t POS_IRS = 2;
constexpr static std::size_t POS_IRE = 3;
constexpr static std::size_t POS_NOM = 9;
constexpr static std::size_t POS_NLE = 6;
constexpr static std::size_t POS_NRS = 7;
constexpr static std::size_t POS_NRE = 8;
constexpr static std::size_t POS_DIR = 4;

namespace lazybastard {

void BlastFileReader::read() {
  threading::WaitGroup wg;
  auto jobFn = [this](threading::Job const *const pJob) { parseLine(pJob); };

  std::string line;
  while (std::getline(m_inputStream, line)) {
    wg.add(1);

    auto job = threading::Job(jobFn, &wg, line);
    m_pThreadPool->addJob(std::move(job));
  }

  wg.wait();
}

void BlastFileReader::parseLine(gsl::not_null<threading::Job const *> const pJob) {
  std::vector<std::string> tokens;

  std::istringstream iss(std::any_cast<std::string>(pJob->getParam(1)), std::ios_base::in);
  std::string token;
  while (std::getline(iss, token, '\t')) {
    tokens.push_back(token);
  }

  if (tokens.size() < std::max({POS_IID, POS_NID, POS_IRS, POS_IRE, POS_NOM, POS_NLE, POS_NRS, POS_NRE, POS_DIR})) {
    throw std::runtime_error("Invalid BLAST file.");
  }

  auto const illuminaRange = std::make_pair(std::stoi(tokens[POS_IRS]), std::stoi(tokens[POS_IRE]) - 1);
  auto const matches = static_cast<std::size_t>(std::stoi(tokens[POS_NOM]));

  auto const nanoporeLength = std::stoi(tokens[POS_NLE]);

  auto addNode = matches >= MINIMUM_MATCHES;
  addNode &= illuminaRange.second - illuminaRange.first + 1 >= static_cast<int>(MINIMUM_MATCHES);

  if (addNode) {
    auto spVertex = std::make_shared<graph::Vertex>(tokens[POS_NID], nanoporeLength);
    m_pGraph->addVertex(std::move(spVertex));

    auto const &nanoporeID = tokens[POS_NID];
    auto const &illuminaID = tokens[POS_IID];

    auto const nanoporeRange = std::make_pair(std::stoi(tokens[POS_NRS]), std::stoi(tokens[POS_NRE]) - 1);
    auto const direction = tokens[POS_DIR] == "+";
    auto const rRatio = static_cast<float>(illuminaRange.second - illuminaRange.first + 1) /
                        static_cast<float>(nanoporeRange.second - nanoporeRange.first + 1);

    auto isPrimary = illuminaRange.second - illuminaRange.first + 1 >= static_cast<int>(TH_LENGTH);
    isPrimary &= matches >= TH_MATCHES;

    auto spVertexMatch = lazybastard::util::make_shared_aggregate<lazybastard::matching::VertexMatch>(
        nanoporeRange, illuminaRange, rRatio, direction, matches, isPrimary);
    m_pMatchMap->addVertexMatch(nanoporeID, illuminaID, std::move(spVertexMatch));
  }

  std::any_cast<threading::WaitGroup *const>(pJob->getParam(0))->done();
}

} // namespace lazybastard

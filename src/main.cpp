#include <algorithm>
#include <any>
#include <cstddef>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <gsl/pointers>
#include <gsl/span>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <lb/BlastFileReader.h>
#include <lb/OutputWriter.h>
#include <lb/Prokrastinator.h>
#include <lb/Util.h>
#include <lb/coroutine/generator.h>
#include <lb/graph/Edge.h>
#include <lb/graph/Graph.h>
#include <lb/graph/Vertex.h>
#include <lb/matching/ID2OverlapMap.h>
#include <lb/matching/MatchMap.h>
#include <lb/threading/Job.h>
#include <lb/threading/ThreadPool.h>
#include <lb/threading/WaitGroup.h>
#include <lb/types/Toggle.h>

#include "Application.h"

//// USING DECLARATIONS ////

using lazybastard::GraphUtil;
using lazybastard::OutputWriter;

using lazybastard::graph::ConsensusDirection;
using lazybastard::graph::Edge;
using lazybastard::graph::EdgeOrder;
using lazybastard::graph::Graph;
using lazybastard::graph::Vertex;

using lazybastard::matching::ID2OverlapMap;
using lazybastard::matching::MatchMap;
using lazybastard::matching::VertexMatch;

using lazybastard::threading::Job;
using lazybastard::threading::ThreadPool;
using lazybastard::threading::WaitGroup;

using lazybastard::util::LTCmp;
using lazybastard::util::make_not_null_and_const;

//// CONSTANTS ////

constexpr static auto BASEWEIGHT_MULTIPLIKATOR = 1.1F;
constexpr static auto MAXWEIGHT_MULTIPLIKATOR = 0.8F;

//// TYPES ////

struct ContainElement { // NOLINT
  std::unordered_map<std::string const *, VertexMatch const *> const matches;
  Vertex const *const nano;
  std::size_t const nanoporeLength;
  std::size_t const score;
  bool const direction;
  bool const isPrimary;
};

//// FUNCTION DECLARATIONS ////

void chainingAndOverlaps(gsl::not_null<Job const *> pJob);
void findContractionEdges(gsl::not_null<Job const *> pJob);
void findContractionTargets(gsl::not_null<Job const *> pJob);
void findDeletableVertices(gsl::not_null<Job const *> pJob);
void contract(gsl::not_null<Job const *> pJob);
void findDeletableEdges(gsl::not_null<Job const *> pJob);
void computeBitweight(gsl::not_null<Job const *> pJob);
void decycle(gsl::not_null<Job const *> pJob);
void assemblePaths(gsl::not_null<Job const *> pJob);

//// FUNCTION DEFINITIONS ////

auto main(int const argc, char const *argv[]) -> int {
  gsl::span<char const *> const args = {argv, static_cast<std::size_t>(argc)};
  Application app(args);

  if (!app.checkIntegrity()) {
    std::cerr << "Paths are pointing to invalid/unusable locations" << std::endl;

    return -1;
  }

  // Keep one free for the main program
  auto threadCount = app.getThreadCount() - 1;
  auto threadPool = ThreadPool(threadCount);

  // Read BLAST file
  auto graph = Graph();
  auto matchMap = MatchMap(&threadPool, &graph);
  if (std::ifstream inputStream{app.getContigsFilePath(), std::ios::binary | std::ios::in}) {
    auto blastReader = lazybastard::BlastFileReader(&threadPool, inputStream, &graph, &matchMap);
    blastReader.read();
  } else {
    std::cerr << "Can't open BLAST file for reading" << std::endl;

    return -1;
  }

  try {
    matchMap.calculateEdges();

    std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

    WaitGroup wg;

    auto chainingJob = [](Job const *const pJob) { chainingAndOverlaps(pJob); };
    auto edges = graph.getEdges();
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);

      auto job = Job(chainingJob, &wg, &matchMap, pEdge, app.getWiggleRoom());
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    std::mutex mutex;
    std::unordered_map<Edge const *, EdgeOrder const *> contractionEdges;
    auto contractionEdgesJob = [](Job const *const pJob) { findContractionEdges(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
      std::for_each(std::begin(pEdge->getEdgeOrders()), std::end(pEdge->getEdgeOrders()), [&](auto const &order) {
        wg.add(1);

        auto job =
            Job(contractionEdgesJob, &wg, &graph, &order, &contractionEdges, std::ref(mutex), app.getWiggleRoom());
        threadPool.addJob(std::move(job));
      });
    });
    wg.wait();

    std::cout << "Number of contraction edges: " << contractionEdges.size() << std::endl;

    std::unordered_map<Vertex const *, Vertex const *> contractionTargets;
    auto contractionTargetsJob = [](Job const *const pJob) { findContractionTargets(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(contractionTargetsJob, &wg, pOrder, &contractionTargets, std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::set<std::string const *const> deletableVertices;
    std::set<Vertex const *const> contractionRoots;
    auto deletableVerticesJob = [](Job const *const pJob) { findDeletableVertices(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(deletableVerticesJob, &wg, pOrder, &deletableVertices, &contractionRoots,
                                             &contractionTargets, std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::unordered_map<Vertex const *, std::unique_ptr<ContainElement const> const> containElements;
    auto contractionJob = [](Job const *const pJob) { contract(pJob); };
    for (auto const &[pEdge, pOrder] : contractionEdges) {
      LB_UNUSED(pEdge);

      wg.add(1);
      auto job = lazybastard::threading::Job(contractionJob, &wg, pOrder, &containElements, &contractionRoots,
                                             &matchMap, std::ref(mutex));
      threadPool.addJob(std::move(job));
    }
    wg.wait();

    std::set<Edge const *const> deletableEdges;
    auto deletableEdgesJob = [](Job const *const pJob) { findDeletableEdges(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);
      auto job = Job(deletableEdgesJob, &wg, pEdge, &deletableEdges, std::ref(mutex));
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    std::cout << "Edges to become deleted: " << deletableEdges.size() << std::endl;

    std::for_each(std::begin(deletableEdges), std::end(deletableEdges),
                  [&](auto const *const pTarget) { graph.deleteEdge(pTarget); });
    deletableEdges.clear();

    std::cout << "Vertices to become deleted: " << deletableVertices.size() << std::endl;

    std::for_each(std::begin(deletableVertices), std::end(deletableVertices),
                  [&](auto const *const pTarget) { graph.deleteVertex(pTarget); });
    deletableVertices.clear();

    std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl;

    edges = graph.getEdges();
    auto computeBitweightJob = [](Job const *const pJob) { computeBitweight(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto *const pEdge) {
      wg.add(1);
      auto job = Job(computeBitweightJob, &wg, pEdge);
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    auto const maxSpanTree = lazybastard::getMaxSpanTree(&graph);

    auto decycleJob = [](Job const *const pJob) { decycle(pJob); };
    std::for_each(std::begin(edges), std::end(edges), [&](auto const *const pEdge) {
      wg.add(1);
      auto job = Job(decycleJob, &wg, &graph, pEdge, maxSpanTree.get(), &deletableEdges, std::ref<std::mutex>(mutex));
      threadPool.addJob(std::move(job));
    });
    wg.wait();

    std::cout << "Edges to become deleted: " << deletableEdges.size() << std::endl;

    std::for_each(std::begin(deletableEdges), std::end(deletableEdges),
                  [&](auto const *const pTarget) { graph.deleteEdge(pTarget); });
    deletableEdges.clear();
    edges = graph.getEdges();

    std::cout << "Order: " << graph.getOrder() << " Size: " << graph.getSize() << std::endl
              << "Starting assembly" << std::endl;

    //// OUTPUT FILES ////
    std::ofstream ofQuery(app.getOutputFilePath() + "/temp.query.fa", std::ios::binary);
    std::ofstream ofPAF(app.getOutputFilePath() + "/temp.align.paf", std::ios::binary);
    std::ofstream ofTarget(app.getOutputFilePath() + "/temp.target.fa", std::ios::binary);

    OutputWriter outputWriter(ofQuery, ofPAF, ofTarget);
    //// /OUTPUT FILES ////

    std::size_t assemblyIdx = 0;
    auto id2OverlapMap = ID2OverlapMap();
    auto connectedComponents = lazybastard::getConnectedComponents(&graph);
    auto assemblyJob = [](Job const *const pJob) { assemblePaths(pJob); };
    for (auto const &connectedComponent : connectedComponents) {

      wg.add(1);
      auto job = Job(assemblyJob, &wg, &graph, &matchMap, &id2OverlapMap, &connectedComponent, assemblyIdx,
                     std::ref(outputWriter));
      threadPool.addJob(std::move(job));

      ++assemblyIdx;
    }
    wg.wait();

    std::cout << "Finished assembly" << std::endl;
  } catch (std::exception const &e) {
    std::cerr << "Exception occurred: " << '\n' << e.what();
  }
  return 0;
}

void chainingAndOverlaps(gsl::not_null<Job const *> const pJob) {
  std::set<gsl::not_null<std::string const *> const, LTCmp<gsl::not_null<std::string const *> const>> plusIDs;
  std::set<gsl::not_null<std::string const *> const, LTCmp<gsl::not_null<std::string const *> const>> minusIDs;

  auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(1)));
  auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(2)));
  auto const *const pEdgeMatches = pMatchMap->getEdgeMatches(pEdge->getID());

  if (!pEdgeMatches) {
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
    return;
  }

  for (auto const &[illuminaID, edgeMatch] : *pEdgeMatches) {
    if (edgeMatch->direction) {
      plusIDs.insert(&illuminaID);
    } else {
      minusIDs.insert(&illuminaID);
    }
  }

  auto const wiggleRoom = std::any_cast<std::size_t>(pJob->getParam(3));
  auto plusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pEdge, plusIDs, true, wiggleRoom);
  auto minusPaths = lazybastard::getMaxPairwisePaths(pMatchMap, pEdge, minusIDs, false, wiggleRoom);

  auto hasPrimary =
      std::any_of(plusPaths.begin(), plusPaths.end(), [](auto const &plusPath) { return std::get<2>(plusPath); });

  if (!hasPrimary) {
    hasPrimary =
        std::any_of(minusPaths.begin(), minusPaths.end(), [](auto const &minusPath) { return std::get<2>(minusPath); });
  }

  if (hasPrimary) {
    plusPaths.erase(
        std::remove_if(plusPaths.begin(), plusPaths.end(), [](auto const &plusPath) { return !std::get<2>(plusPath); }),
        plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](auto const &minusPath) { return !std::get<2>(minusPath); }),
                     minusPaths.end());
  }

  auto hasMulti = std::any_of(plusPaths.begin(), plusPaths.end(),
                              [](auto const &plusPath) { return std::get<0>(plusPath).size() > 1; });

  if (!hasMulti) {
    hasMulti = std::any_of(minusPaths.begin(), minusPaths.end(),
                           [](auto const &minusPath) { return std::get<0>(minusPath).size() > 1; });
  }

  if (hasMulti) {
    plusPaths.erase(std::remove_if(plusPaths.begin(), plusPaths.end(),
                                   [](auto const &plusPath) { return std::get<0>(plusPath).size() <= 1; }),
                    plusPaths.end());
    minusPaths.erase(std::remove_if(minusPaths.begin(), minusPaths.end(),
                                    [](auto const &minusPath) { return std::get<0>(minusPath).size() <= 1; }),
                     minusPaths.end());
  }

  auto const combinedSize = plusPaths.size() + minusPaths.size();
  if (combinedSize > 1) {
    pEdge->setShadow(true);
  } else if (combinedSize == 1) {
    auto const &path = !minusPaths.empty() ? minusPaths[0] : plusPaths[0];
    pEdge->setShadow(std::get<2>(path));
  }

  std::for_each(std::begin(plusPaths), std::end(plusPaths), [=](auto &plusPath) {
    auto const overlap = lazybastard::computeOverlap(pMatchMap, std::move(std::get<0>(plusPath)), pEdge, false,
                                                     std::get<1>(plusPath), std::get<2>(plusPath));
    if (overlap.has_value()) {
      pEdge->appendOrder(std::move(overlap.value()));
    }
  });

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findContractionEdges(gsl::not_null<Job const *> const pJob) {
  auto const pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
  auto const pOrder = make_not_null_and_const(std::any_cast<EdgeOrder const *>(pJob->getParam(2)));
  if (pOrder->isContained && pOrder->isPrimary) {
    auto isSane = true;
    auto const *const pNeighbors = pGraph->getNeighbors(pOrder->startVertex->getID());
    auto const edges = std::map<std::string, Edge *>(std::begin(*pNeighbors), std::end(*pNeighbors));
    for (auto const &[targetID, pEdge] : edges) {
      if (targetID == pOrder->endVertex->getID() || pEdge->isShadow()) {
        continue;
      }
      isSane &= pGraph->hasEdge(std::make_pair(make_not_null_and_const(&pOrder->endVertex->getID()), &targetID));
      if (!isSane) {
        break;
      }

      auto const *const pTargetVertex =
          pEdge->getVertices().first != pOrder->startVertex ? pEdge->getVertices().first : pEdge->getVertices().second;

      auto const wiggleRoom = std::any_cast<std::size_t>(pJob->getParam(5));

      isSane &=
          lazybastard::sanityCheck(pGraph, pOrder->startVertex, pOrder->endVertex, pTargetVertex, pOrder, wiggleRoom);

      if (!isSane) {
        break;
      }

      std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(4)).get());
      auto const pContractionEdges =
          gsl::make_not_null(std::any_cast<std::unordered_map<Edge const *, EdgeOrder const *> *>(pJob->getParam(3)));
      pContractionEdges->insert({pEdge, pOrder});
    }
  }
  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findContractionTargets(gsl::not_null<Job const *> const pJob) {
  auto const pOrder = gsl::make_not_null(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pContractionTargets =
      gsl::make_not_null(std::any_cast<std::unordered_map<Vertex const *, Vertex const *> *>(pJob->getParam(2)));
  {
    std::lock_guard<std::mutex> guard(std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get());

    auto const iter = pContractionTargets->find(pOrder->endVertex);
    auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->endVertex : iter->second;

    if (pContractionTargets->find(pOrder->startVertex) == pContractionTargets->end()) {
      pContractionTargets->insert({pOrder->startVertex, contractTo});
    }
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableVertices(gsl::not_null<const Job *> pJob) {
  auto const pOrder = gsl::make_not_null(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pDeletableVertices =
      gsl::make_not_null(std::any_cast<std::set<std::string const *const> *>(pJob->getParam(2)));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
    pDeletableVertices->insert(&pOrder->startVertex->getID());
  }

  auto const pContractionTargets =
      gsl::make_not_null(std::any_cast<std::unordered_map<Vertex const *, Vertex const *> *>(pJob->getParam(4)));
  auto const pContractionRoots = gsl::make_not_null(std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(3)));
  {
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT

    auto const iter = pContractionTargets->find(pOrder->startVertex);
    auto const *const contractTo = iter == pContractionTargets->end() ? pOrder->startVertex : iter->second;

    pContractionRoots->insert(contractTo);
    pContractionRoots->erase(pOrder->startVertex);
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void contract(gsl::not_null<Job const *> const pJob) {
  auto const pOrder = make_not_null_and_const(std::any_cast<EdgeOrder const *>(pJob->getParam(1)));
  auto const pContractionRoots = gsl::make_not_null(std::any_cast<std::set<Vertex const *const> *>(pJob->getParam(3)));
  auto const iter = pContractionRoots->find(pOrder->endVertex);

  if (iter == pContractionRoots->end()) {
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
    return;
  }

  auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(4)));
  std::unordered_map<std::string const *, VertexMatch const *> matches;
  std::for_each(std::begin(pOrder->ids), std::end(pOrder->ids), [&](auto const id) {
    auto const *pVertexMatches = pMatchMap->getVertexMatch(pOrder->startVertex->getID(), *id);
    if (pVertexMatches != nullptr) {
      matches.insert({id, pVertexMatches});
    }
  });

  std::lock_guard<std::mutex> guard(
      std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
  auto const pContainElements = gsl::make_not_null(
      std::any_cast<std::unordered_map<Vertex const *, std::unique_ptr<ContainElement const> const> *>(
          pJob->getParam(2)));
  pContainElements->insert(std::make_pair(
      pOrder->endVertex, std::make_unique<ContainElement>(ContainElement{
                             std::move(matches), pOrder->startVertex, pOrder->startVertex->getNanoporeLength(),
                             pOrder->score, pOrder->direction, pOrder->isPrimary})));

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void findDeletableEdges(gsl::not_null<Job const *> const pJob) {
  auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
  std::vector<EdgeOrder> filteredOrders;
  std::remove_copy_if(std::begin(pEdge->getEdgeOrders()), std::end(pEdge->getEdgeOrders()),
                      std::back_inserter(filteredOrders), [](auto const &order) { return order.isContained; });

  if (filteredOrders.empty()) {
    auto const pDeletableEdges = gsl::make_not_null(std::any_cast<std::set<Edge const *const> *>(pJob->getParam(2)));
    std::lock_guard<std::mutex> guard(
        std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(3)).get()); // NOLINT
    pDeletableEdges->insert(pEdge.get());
  }

  pEdge->replaceOrders(std::move(filteredOrders));
  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void computeBitweight(gsl::not_null<Job const *> const pJob) {
  auto const pEdge = gsl::make_not_null(std::any_cast<Edge *>(pJob->getParam(1)));
  auto const &orders = pEdge->getEdgeOrders();

  if (orders.empty()) {
    std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
    return;
  }

  if (pEdge->isShadow()) {
    auto const initial = orders[0].direction;
    auto const findOther = std::find_if(orders.begin(), orders.end(),
                                        [&](auto const &order) { return order.direction != initial; }) != orders.end();
    if (!findOther) {
      pEdge->setConsensusDirection(orders[0].direction);
    }
  } else {
    pEdge->setWeight(orders[0].score);
    pEdge->setConsensusDirection(orders[0].direction);
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void decycle(gsl::not_null<Job const *> const pJob) {
  auto const pEdge = gsl::make_not_null(std::any_cast<Edge const *>(pJob->getParam(2)));
  auto const pMaxSpanTree = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(3)));
  if (pEdge->getConsensusDirection() != ConsensusDirection::Enum::e_NONE &&
      !pMaxSpanTree->hasEdge(std::make_pair(make_not_null_and_const(&pEdge->getVertices().first->getID()),
                                            make_not_null_and_const(&pEdge->getVertices().second->getID())))) {
    auto const pGraph = make_not_null_and_const(std::any_cast<Graph *>(pJob->getParam(1)));
    auto const shortestPath = GraphUtil::getShortestPath(pGraph, pEdge->getVertices());
    lazybastard::Toggle direction = pEdge->getConsensusDirection();
    auto const baseWeight = static_cast<float>(pEdge->getWeight());
    std::vector<float> weights;

    for (auto it = shortestPath.begin(); it != std::prev(shortestPath.end()); ++it) {
      auto const *const pPathEdge = pGraph->getEdge(std::make_pair(&(*it)->getID(), &(*std::next(it))->getID()));
      if (pPathEdge != nullptr) {
        direction = direction && pPathEdge->getConsensusDirection();
        weights.push_back(static_cast<float>(pEdge->getWeight()));
      }
    }

    if (!direction && !weights.empty()) {
      auto const iterMinWeight = std::min_element(weights.begin(), weights.end());
      auto const iterMaxWeight = std::max_element(weights.begin(), weights.end());

      if (*iterMaxWeight < baseWeight || (baseWeight * BASEWEIGHT_MULTIPLIKATOR >= *iterMinWeight &&
                                          *iterMinWeight < *iterMaxWeight * MAXWEIGHT_MULTIPLIKATOR)) {
        auto const minWeightIdx = std::distance(weights.begin(), iterMinWeight);
        auto const *const pDeletableEdge =
            pGraph->getEdge(std::make_pair(&(*std::next(shortestPath.begin(), minWeightIdx))->getID(),
                                           &(*std::next(shortestPath.begin(), minWeightIdx + 1))->getID()));
        auto const pDeletableEdges =
            gsl::make_not_null(std::any_cast<std::set<Edge const *const> *>(pJob->getParam(4)));
        std::lock_guard<std::mutex> guard(
            std::any_cast<std::reference_wrapper<std::mutex>>(pJob->getParam(5)).get()); // NOLINT
        pDeletableEdges->insert(pDeletableEdge);
      }
    }
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}

void assemblePaths(gsl::not_null<Job const *> const pJob) {
  auto const pGraph = gsl::make_not_null(std::any_cast<Graph *>(pJob->getParam(1)));
  auto const pConnectedComponent =
      make_not_null_and_const(std::any_cast<std::vector<lazybastard::graph::Vertex *> const *>(pJob->getParam(4)));

  auto const pSubGraph = pGraph->getSubgraph(*pConnectedComponent);
  auto const subGraphVertices = pSubGraph->getVertices();
  auto const maxNPLVertexIter = std::max_element(
      std::begin(subGraphVertices), std::end(subGraphVertices),
      [](Vertex const *v1, Vertex const *v2) { return v1->getNanoporeLength() < v2->getNanoporeLength(); });
  auto *const pMaxNPLVertex = maxNPLVertexIter == pSubGraph->getVertices().end() ? nullptr : *maxNPLVertexIter;
  auto const pMatchMap = make_not_null_and_const(std::any_cast<MatchMap *>(pJob->getParam(2)));

  if (pMaxNPLVertex != nullptr) {
    auto const pDiGraph = lazybastard::getDirectionGraph(pGraph, pSubGraph.get(), pMaxNPLVertex);
    auto const paths = lazybastard::linearizeGraph(pDiGraph.get());

    auto const pID2OverlapMap = gsl::make_not_null(std::any_cast<ID2OverlapMap *>(pJob->getParam(3)));
    auto const assemblyIdx = std::any_cast<std::size_t>(pJob->getParam(5));
    std::for_each(paths.begin(), paths.end(), [&](auto const &path) {
      lazybastard::assemblePath(pGraph, pMatchMap, pID2OverlapMap, &path, pDiGraph.get(), assemblyIdx,
                                std::any_cast<std::reference_wrapper<OutputWriter>>(pJob->getParam(6)).get());
    });
  }

  std::any_cast<WaitGroup *>(pJob->getParam(0))->done();
}
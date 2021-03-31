#include "Prokrastinator.h"

#include <utility>

#include "Util.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "matching/MatchMap.h"

namespace {

std::pair<double, double> computeOverhangs(gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
                                           gsl::not_null<lazybastard::graph::Vertex const *> const pVertex,
                                           gsl::not_null<lazybastard::graph::Edge const *> const pEdge,
                                           std::string const &illuminaID) {
  auto const pVertexMatch = gsl::make_not_null(pMatches->getVertexMatch(pVertex->getID(), illuminaID));
  auto const pEdgeMatch = gsl::make_not_null(pMatches->getEdgeMatch(pEdge->getID(), illuminaID));

  auto nanoCorrectionLeft = static_cast<double>(pEdgeMatch->overlap.first - pVertexMatch->illuminaRange.first) //
                            / pVertexMatch->rRatio;
  auto nanoCorrectionRight = static_cast<double>(pVertexMatch->illuminaRange.second - pEdgeMatch->overlap.second) //
                             / pVertexMatch->rRatio;

  lazybastard::util::swap_if(nanoCorrectionLeft, nanoCorrectionRight, !pVertexMatch->direction);

  auto const overhangLeft = static_cast<double>(pVertexMatch->nanoporeRange.first) + nanoCorrectionLeft;
  auto const nanoporeLength = static_cast<int>(pVertex->getNanoporeLength());
  auto const overhangRight = static_cast<double>(nanoporeLength - pVertexMatch->nanoporeRange.second) //
                             + nanoCorrectionRight;

  return std::make_pair(overhangLeft, overhangRight);
}

} // unnamed namespace

std::optional<lazybastard::graph::EdgeOrder>
lazybastard::computeOverlap(gsl::not_null<lazybastard::matching::MatchMap const *> const pMatches,
                            std::deque<std::string> &ids, gsl::not_null<lazybastard::graph::Edge const *> const pEdge,
                            bool direction, std::size_t score, bool isPrimary) {
  auto const &firstID = ids.front();
  auto const &lastID = ids.back();

  auto const vertices = pEdge->getVertices();

  auto const overhangsFirstIDFirstVertex = computeOverhangs(pMatches, vertices.first, pEdge, firstID);
  auto const overhangsLastIDFirstVertex = computeOverhangs(pMatches, vertices.first, pEdge, lastID);
  auto const overhangsFirstIDSecondVertex = computeOverhangs(pMatches, vertices.second, pEdge, firstID);
  auto const overhangsLastIDSecondVertex = computeOverhangs(pMatches, vertices.second, pEdge, lastID);

  auto const leftOverhangFirstVertex = overhangsFirstIDFirstVertex.first;
  auto const rightOverhangFirstVertex = overhangsLastIDFirstVertex.second;

  auto leftOverhangSecondVertex = overhangsFirstIDSecondVertex.first;
  auto rightOverhangSecondVertex = overhangsLastIDSecondVertex.second;

  if (!direction) {
    leftOverhangSecondVertex = overhangsFirstIDSecondVertex.second;
    rightOverhangSecondVertex = overhangsLastIDSecondVertex.first;
  }

  std::optional<graph::EdgeOrder> result;
  if (leftOverhangFirstVertex <= leftOverhangSecondVertex && rightOverhangFirstVertex <= rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.first, vertices.second, leftOverhangSecondVertex - leftOverhangFirstVertex,
                                    rightOverhangSecondVertex - rightOverhangFirstVertex, true, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex >= leftOverhangSecondVertex &&
             rightOverhangFirstVertex >= rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.second, vertices.first, leftOverhangFirstVertex - leftOverhangSecondVertex,
                                    rightOverhangFirstVertex - rightOverhangSecondVertex, true, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex > leftOverhangSecondVertex &&
             rightOverhangFirstVertex < rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.first, vertices.second, leftOverhangFirstVertex - leftOverhangSecondVertex,
                                    rightOverhangSecondVertex - rightOverhangFirstVertex, false, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  } else if (leftOverhangFirstVertex < leftOverhangSecondVertex &&
             rightOverhangFirstVertex > rightOverhangSecondVertex) {
    result.emplace(graph::EdgeOrder{vertices.second, vertices.first, leftOverhangSecondVertex - leftOverhangFirstVertex,
                                    rightOverhangFirstVertex - rightOverhangSecondVertex, false, vertices.first, score,
                                    std::move(ids), direction, isPrimary});
  }

  return result;
}
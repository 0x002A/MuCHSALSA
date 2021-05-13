#include <gsl/pointers>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "Util.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

TEST(GraphTest, BasicTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::make_pair(spVertex1.get(), spVertex2.get())), "1,2");

  auto vertexIDs = std::make_pair(spVertex1.get(), spVertex2.get());
  auto graph     = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertexIDs);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::move(vertexIDs)), "1,2");

  auto expected    = std::make_pair(spVertex1.get(), spVertex2.get());
  auto expectedToo = std::make_pair(spVertex2.get(), spVertex1.get());
  ASSERT_TRUE(graph.hasEdge(expected));
  ASSERT_TRUE(graph.hasEdge(expectedToo));
}

TEST(GraphTest, EdgeCreationTest) {
  auto const spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto const spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);

  auto const edge = lazybastard::graph::Edge(std::make_pair(spVertex1, spVertex2));

  auto const vertices = edge.getVertices();
  ASSERT_EQ(vertices.first->getId(), "2");
  ASSERT_EQ(vertices.second->getId(), "1");
}

TEST(GraphTest, InverseOrderTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::make_pair(spVertex1.get(), spVertex2.get())), "2,1");

  auto vertexIDs = std::make_pair(spVertex1.get(), spVertex2.get());
  auto graph     = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertexIDs);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::move(vertexIDs)), "2,1");

  auto expected    = std::make_pair(spVertex1.get(), spVertex2.get());
  auto expectedToo = std::make_pair(spVertex2.get(), spVertex1.get());
  ASSERT_TRUE(graph.hasEdge(expected));
  ASSERT_TRUE(graph.hasEdge(expectedToo));
}

TEST(GraphTest, DirectedGraphTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::make_pair(spVertex1.get(), spVertex2.get())), "2,1");

  auto vertices = std::make_pair(spVertex1.get(), spVertex2.get());
  auto graph    = lazybastard::graph::DiGraph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertices);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeId(std::move(vertices)), "2,1");

  auto expected = std::make_pair(gsl::make_not_null(spVertex1.get()), gsl::make_not_null(spVertex2.get()));
  ASSERT_TRUE(graph.hasEdge(expected));

  auto notExpected = std::make_pair(gsl::make_not_null(spVertex2.get()), gsl::make_not_null(spVertex1.get()));
  ASSERT_FALSE(graph.hasEdge(notExpected));
}

TEST(GraphTest, EdgeDeletionTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);

  auto firstEdge  = std::make_pair(spVertex1.get(), spVertex2.get());
  auto secondEdge = std::make_pair(spVertex1.get(), spVertex3.get());
  auto thirdEdge  = std::make_pair(spVertex2.get(), spVertex3.get());

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 3);

  graph.addEdge(firstEdge);
  graph.addEdge(secondEdge);
  graph.addEdge(thirdEdge);
  ASSERT_EQ(graph.getSize(), 3);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(std::move(spVertex1));
  diGraph.addVertex(std::move(spVertex2));
  diGraph.addVertex(std::move(spVertex3));
  ASSERT_EQ(diGraph.getOrder(), 3);

  diGraph.addEdge(firstEdge);
  diGraph.addEdge(secondEdge);
  diGraph.addEdge(thirdEdge);
  ASSERT_EQ(diGraph.getSize(), 3);

  auto *pSecondEdge = graph.getEdge(std::make_pair(secondEdge.first, secondEdge.second));
  ASSERT_NE(pSecondEdge, nullptr);
  graph.deleteEdge(pSecondEdge, nullptr);
  ASSERT_EQ(graph.getOrder(), 3);
  ASSERT_EQ(graph.getSize(), 2);

  ASSERT_TRUE(graph.hasEdge(std::make_pair(firstEdge.first, firstEdge.second)));
  ASSERT_TRUE(graph.hasEdge(std::make_pair(thirdEdge.first, thirdEdge.second)));

  pSecondEdge = diGraph.getEdge(std::make_pair(secondEdge.first, secondEdge.second));
  ASSERT_NE(pSecondEdge, nullptr);
  diGraph.deleteEdge(pSecondEdge, nullptr);
  ASSERT_EQ(diGraph.getOrder(), 3);
  ASSERT_EQ(diGraph.getSize(), 2);

  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(firstEdge.first, firstEdge.second)));
  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(thirdEdge.first, thirdEdge.second)));
}

TEST(GraphTest, VertexDeletionTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 3);

  graph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  graph.addEdge(std::make_pair(spVertex2.get(), spVertex3.get()));
  ASSERT_EQ(graph.getSize(), 3);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  ASSERT_EQ(diGraph.getOrder(), 3);

  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  diGraph.addEdge(std::make_pair(spVertex2.get(), spVertex3.get()));
  ASSERT_EQ(diGraph.getSize(), 3);

  graph.deleteVertex(spVertex2.get(), nullptr);
  ASSERT_EQ(graph.getOrder(), 2);
  ASSERT_EQ(graph.getSize(), 1);

  ASSERT_TRUE(graph.hasEdge(std::make_pair(spVertex1.get(), spVertex3.get())));

  diGraph.deleteVertex(spVertex2.get(), nullptr);
  ASSERT_EQ(diGraph.getOrder(), 2);
  ASSERT_EQ(diGraph.getSize(), 1);

  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(spVertex1.get(), spVertex3.get())));
}

TEST(GraphTest, NeighboorTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>("4", 0);

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  graph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 4);

  graph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  ASSERT_EQ(graph.getSize(), 4);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  diGraph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  ASSERT_EQ(diGraph.getSize(), 4);

  auto graphNeighbors = graph.getNeighbors(spVertex3.get());
  ASSERT_EQ(graphNeighbors.size(), 3);
  ASSERT_TRUE(graphNeighbors.contains("1"));
  ASSERT_TRUE(graphNeighbors.contains("2"));
  ASSERT_TRUE(graphNeighbors.contains("4"));

  graphNeighbors = diGraph.getPredecessors(spVertex3.get());
  ASSERT_EQ(graphNeighbors.size(), 1);
  ASSERT_TRUE(graphNeighbors.contains("1"));

  auto const diGraphNeighbors = diGraph.getSuccessors(spVertex3.get());
  ASSERT_EQ(diGraphNeighbors.size(), 2);
  ASSERT_TRUE(diGraphNeighbors.contains("2"));
  ASSERT_TRUE(diGraphNeighbors.contains("4"));
}

TEST(GraphTest, SubgraphTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>("4", 0);

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  graph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 4);

  graph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  ASSERT_EQ(graph.getSize(), 4);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  diGraph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  ASSERT_EQ(diGraph.getSize(), 4);

  auto       v1 = spVertex1.get(), v2 = spVertex3.get();
  auto const subgraphVertices = std::vector<lazybastard::graph::Vertex *>({v1, v2});
  auto const subgraph         = graph.getSubgraph(subgraphVertices);
  ASSERT_EQ(subgraph->getSize(), 1);
  ASSERT_TRUE(subgraph->hasVertex(v1->getId()));
  ASSERT_TRUE(subgraph->hasVertex(v2->getId()));

  ASSERT_TRUE(subgraph->hasEdge(std::make_pair(v1, v2)));
  ASSERT_TRUE(subgraph->hasEdge(std::make_pair(v2, v1)));

  auto const subDigraph = diGraph.getSubgraph(subgraphVertices);
  ASSERT_EQ(subgraph->getSize(), 1);
  ASSERT_TRUE(subDigraph->hasVertex(v1->getId()));
  ASSERT_TRUE(subDigraph->hasVertex(v2->getId()));

  ASSERT_TRUE(subDigraph->hasEdge(std::make_pair(v1, v2)));

  auto const &inDegrees = subDigraph->getInDegrees();
  ASSERT_EQ(inDegrees.size(), 2);
  ASSERT_EQ(inDegrees.at(v1), 0);
  ASSERT_EQ(inDegrees.at(v2), 1);

  auto const &outDegrees = subDigraph->getOutDegrees();
  ASSERT_EQ(outDegrees.size(), 2);
  ASSERT_EQ(outDegrees.at(v1), 1);
  ASSERT_EQ(outDegrees.at(v2), 0);
}

TEST(GraphTest, ShortestPathTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>("4", 0);

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  graph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(graph.getOrder(), 4);

  graph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex2.get(), spVertex3.get()));
  graph.addEdge(std::make_pair(spVertex2.get(), spVertex4.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  graph.addEdge(std::make_pair(spVertex4.get(), spVertex1.get()));
  ASSERT_EQ(graph.getSize(), 5);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  diGraph.addVertex(spVertex4->getSharedPtr());
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  diGraph.addEdge(std::make_pair(spVertex2.get(), spVertex3.get()));
  diGraph.addEdge(std::make_pair(spVertex2.get(), spVertex4.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  diGraph.addEdge(std::make_pair(spVertex4.get(), spVertex1.get()));
  ASSERT_EQ(diGraph.getSize(), 5);

  std::string v1("1"), v2("4");
  auto        shortestPathVertices = std::make_pair(lazybastard::util::make_not_null_and_const(graph.getVertex(v1)),
                                             lazybastard::util::make_not_null_and_const(graph.getVertex(v2)));
  auto        shortestPath =
      lazybastard::GraphUtil::getShortestPath(lazybastard::util::make_not_null_and_const(&graph), shortestPathVertices);
  ASSERT_EQ(shortestPath.size(), 2);
  ASSERT_EQ(shortestPath[0]->getId(), "1");
  ASSERT_EQ(shortestPath[1]->getId(), "4");

  shortestPathVertices = std::make_pair(lazybastard::util::make_not_null_and_const(diGraph.getVertex(v1)),
                                        lazybastard::util::make_not_null_and_const(diGraph.getVertex(v2)));
  shortestPath         = lazybastard::GraphUtil::getShortestPath(lazybastard::util::make_not_null_and_const(&diGraph),
                                                         shortestPathVertices);
  ASSERT_EQ(shortestPath.size(), 3);
  ASSERT_EQ(shortestPath[0]->getId(), "1");
  ASSERT_EQ(shortestPath[1]->getId(), "2");
  ASSERT_EQ(shortestPath[2]->getId(), "4");
}

TEST(GraphTest, DegreeTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);

  auto firstEdge  = std::make_pair(spVertex1.get(), spVertex2.get());
  auto secondEdge = std::make_pair(spVertex1.get(), spVertex3.get());
  auto thirdEdge  = std::make_pair(spVertex2.get(), spVertex3.get());

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  ASSERT_EQ(diGraph.getOrder(), 3);

  diGraph.addEdge(firstEdge);
  diGraph.addEdge(secondEdge);
  diGraph.addEdge(thirdEdge);
  diGraph.addEdge(thirdEdge); // Test double insertion
  ASSERT_EQ(diGraph.getSize(), 3);

  auto inDegrees = diGraph.getInDegrees();
  ASSERT_EQ(inDegrees.size(), 3);
  ASSERT_EQ(inDegrees.at(spVertex1.get()), 0);
  ASSERT_EQ(inDegrees.at(spVertex2.get()), 1);
  ASSERT_EQ(inDegrees.at(spVertex3.get()), 2);

  auto outDegrees = diGraph.getOutDegrees();
  ASSERT_EQ(outDegrees.size(), 3);
  ASSERT_EQ(outDegrees.at(spVertex1.get()), 2);
  ASSERT_EQ(outDegrees.at(spVertex2.get()), 1);
  ASSERT_EQ(outDegrees.at(spVertex3.get()), 0);

  diGraph.deleteVertex(spVertex2.get(), nullptr);

  inDegrees = diGraph.getInDegrees();
  ASSERT_EQ(inDegrees.size(), 2);
  ASSERT_EQ(inDegrees.at(spVertex1.get()), 0);
  ASSERT_EQ(inDegrees.at(spVertex3.get()), 1);

  outDegrees = diGraph.getOutDegrees();
  ASSERT_EQ(outDegrees.size(), 2);
  ASSERT_EQ(outDegrees.at(spVertex1.get()), 1);
  ASSERT_EQ(outDegrees.at(spVertex3.get()), 0);

  auto *pSecondEdge = diGraph.getEdge(std::make_pair(secondEdge.first, secondEdge.second));
  ASSERT_NE(pSecondEdge, nullptr);
  diGraph.deleteEdge(pSecondEdge, nullptr);

  inDegrees = diGraph.getInDegrees();
  ASSERT_EQ(inDegrees.size(), 2);
  ASSERT_EQ(inDegrees.at(spVertex1.get()), 0);
  ASSERT_EQ(inDegrees.at(spVertex3.get()), 0);

  outDegrees = diGraph.getOutDegrees();
  ASSERT_EQ(outDegrees.size(), 2);
  ASSERT_EQ(outDegrees.at(spVertex1.get()), 0);
  ASSERT_EQ(outDegrees.at(spVertex3.get()), 0);
}

TEST(GraphTest, TopologicalSortTest) {
  auto diGraph = lazybastard::graph::DiGraph();

  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>("4", 0);
  auto spVertex5 = std::make_shared<lazybastard::graph::Vertex>("5", 0);

  diGraph.addVertex(spVertex1->getSharedPtr());
  diGraph.addVertex(spVertex2->getSharedPtr());
  diGraph.addVertex(spVertex3->getSharedPtr());
  diGraph.addVertex(spVertex4->getSharedPtr());
  diGraph.addVertex(spVertex5->getSharedPtr());

  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex4.get()));
  diGraph.addEdge(std::make_pair(spVertex1.get(), spVertex3.get()));
  diGraph.addEdge(std::make_pair(spVertex2.get(), spVertex1.get()));
  diGraph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  diGraph.addEdge(std::make_pair(spVertex5.get(), spVertex2.get()));

  auto const sortedGraph = diGraph.sortTopologically();

  ASSERT_EQ(sortedGraph.size(), 5);

  ASSERT_EQ(sortedGraph[0], spVertex5.get());
  ASSERT_EQ(sortedGraph[1], spVertex2.get());
  ASSERT_EQ(sortedGraph[2], spVertex1.get());
  ASSERT_EQ(sortedGraph[3], spVertex3.get());
  ASSERT_EQ(sortedGraph[4], spVertex4.get());
}
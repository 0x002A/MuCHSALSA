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
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::make_pair(spVertex1->getID(), spVertex2->getID())), "1,2");

  auto vertexIDs = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto graph = lazybastard::graph::Graph();
  graph.addVertex(std::move(spVertex1));
  graph.addVertex(std::move(spVertex2));
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertexIDs);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::move(vertexIDs)), "1,2");

  auto const v1 = std::string("1");
  auto const v2 = std::string("2");
  auto expected = std::make_pair(gsl::make_not_null(&v1), gsl::make_not_null(&v2));
  auto expectedToo = std::make_pair(gsl::make_not_null(&v2), gsl::make_not_null(&v1));
  ASSERT_TRUE(graph.hasEdge(expected));
  ASSERT_TRUE(graph.hasEdge(expectedToo));
}

TEST(GraphTest, EdgeCreationTest) {
  auto const spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto const spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);

  auto const edge = lazybastard::graph::Edge(std::make_pair(spVertex1, spVertex2));

  auto const vertices = edge.getVertices();
  ASSERT_EQ(vertices.first->getID(), "2");
  ASSERT_EQ(vertices.second->getID(), "1");
}

TEST(GraphTest, InverseOrderTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::make_pair(spVertex1->getID(), spVertex2->getID())), "2,1");

  auto vertexIDs = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto graph = lazybastard::graph::Graph();
  graph.addVertex(std::move(spVertex1));
  graph.addVertex(std::move(spVertex2));
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertexIDs);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::move(vertexIDs)), "1,2");

  auto const v1 = std::string("1");
  auto const v2 = std::string("2");
  auto expected = std::make_pair(gsl::make_not_null(&v1), gsl::make_not_null(&v2));
  auto expectedToo = std::make_pair(gsl::make_not_null(&v2), gsl::make_not_null(&v1));
  ASSERT_TRUE(graph.hasEdge(expected));
  ASSERT_TRUE(graph.hasEdge(expectedToo));
}

TEST(GraphTest, DirectedGraphTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::make_pair(spVertex1->getID(), spVertex2->getID())), "2,1");

  auto vertexIDs = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto graph = lazybastard::graph::DiGraph();
  graph.addVertex(std::move(spVertex1));
  graph.addVertex(std::move(spVertex2));
  ASSERT_EQ(graph.getOrder(), 2);

  graph.addEdge(vertexIDs);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::move(vertexIDs)), "2,1");

  auto const v1 = std::string("2");
  auto const v2 = std::string("1");
  auto expected = std::make_pair(gsl::make_not_null(&v1), gsl::make_not_null(&v2));
  ASSERT_TRUE(graph.hasEdge(expected));

  auto notExpected = std::make_pair(gsl::make_not_null(&v2), gsl::make_not_null(&v1));
  ASSERT_FALSE(graph.hasEdge(notExpected));
}

TEST(GraphTest, EdgeDeletionTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);

  auto firstEdge = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto secondEdge = std::make_pair(spVertex1->getID(), spVertex3->getID());
  auto thirdEdge = std::make_pair(spVertex2->getID(), spVertex3->getID());

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

  auto pSecondEdge = graph.getEdge(std::make_pair(&secondEdge.first, &secondEdge.second));
  ASSERT_NE(pSecondEdge, nullptr);
  graph.deleteEdge(pSecondEdge);
  ASSERT_EQ(graph.getOrder(), 3);
  ASSERT_EQ(graph.getSize(), 2);

  ASSERT_TRUE(graph.hasEdge(std::make_pair(&firstEdge.first, &firstEdge.second)));
  ASSERT_TRUE(graph.hasEdge(std::make_pair(&thirdEdge.first, &thirdEdge.second)));

  pSecondEdge = diGraph.getEdge(std::make_pair(&secondEdge.first, &secondEdge.second));
  ASSERT_NE(pSecondEdge, nullptr);
  diGraph.deleteEdge(pSecondEdge);
  ASSERT_EQ(diGraph.getOrder(), 3);
  ASSERT_EQ(diGraph.getSize(), 2);

  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(&firstEdge.first, &firstEdge.second)));
  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(&thirdEdge.first, &thirdEdge.second)));
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

  graph.addEdge(std::make_pair("1", "2"));
  graph.addEdge(std::make_pair("1", "3"));
  graph.addEdge(std::make_pair("2", "3"));
  ASSERT_EQ(graph.getSize(), 3);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(std::move(spVertex1));
  diGraph.addVertex(std::move(spVertex2));
  diGraph.addVertex(std::move(spVertex3));
  ASSERT_EQ(diGraph.getOrder(), 3);

  diGraph.addEdge(std::make_pair("1", "2"));
  diGraph.addEdge(std::make_pair("1", "3"));
  diGraph.addEdge(std::make_pair("2", "3"));
  ASSERT_EQ(diGraph.getSize(), 3);

  std::string vertex2Delete("2"), firstVertex("1"), thirdVertex("3");
  auto const pVertex2Delete = gsl::make_not_null(&vertex2Delete);

  graph.deleteVertex(pVertex2Delete);
  ASSERT_EQ(graph.getOrder(), 2);
  ASSERT_EQ(graph.getSize(), 1);

  ASSERT_TRUE(graph.hasEdge(std::make_pair(&firstVertex, &thirdVertex)));

  diGraph.deleteVertex(pVertex2Delete);
  ASSERT_EQ(diGraph.getOrder(), 2);
  ASSERT_EQ(diGraph.getSize(), 1);

  ASSERT_TRUE(diGraph.hasEdge(std::make_pair(&firstVertex, &thirdVertex)));
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

  graph.addEdge(std::make_pair("1", "2"));
  graph.addEdge(std::make_pair("1", "3"));
  graph.addEdge(std::make_pair("3", "2"));
  graph.addEdge(std::make_pair("3", "4"));
  ASSERT_EQ(graph.getSize(), 4);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(std::move(spVertex1));
  diGraph.addVertex(std::move(spVertex2));
  diGraph.addVertex(std::move(spVertex3));
  diGraph.addVertex(std::move(spVertex4));
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair("1", "2"));
  diGraph.addEdge(std::make_pair("1", "3"));
  diGraph.addEdge(std::make_pair("3", "2"));
  diGraph.addEdge(std::make_pair("3", "4"));
  ASSERT_EQ(diGraph.getSize(), 4);

  auto graphNeighboors = graph.getNeighbors("3");
  ASSERT_EQ(graphNeighboors.size(), 3);
  ASSERT_TRUE(graphNeighboors.contains("1"));
  ASSERT_TRUE(graphNeighboors.contains("2"));
  ASSERT_TRUE(graphNeighboors.contains("4"));

  graphNeighboors = diGraph.getPredecessors("3");
  ASSERT_EQ(graphNeighboors.size(), 1);
  ASSERT_TRUE(graphNeighboors.contains("1"));

  graphNeighboors = diGraph.getSuccessors("3");
  ASSERT_EQ(graphNeighboors.size(), 2);
  ASSERT_TRUE(graphNeighboors.contains("2"));
  ASSERT_TRUE(graphNeighboors.contains("4"));
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

  graph.addEdge(std::make_pair("1", "2"));
  graph.addEdge(std::make_pair("1", "3"));
  graph.addEdge(std::make_pair("3", "2"));
  graph.addEdge(std::make_pair("3", "4"));
  ASSERT_EQ(graph.getSize(), 4);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(std::move(spVertex1));
  diGraph.addVertex(std::move(spVertex2));
  diGraph.addVertex(std::move(spVertex3));
  diGraph.addVertex(std::move(spVertex4));
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair("1", "2"));
  diGraph.addEdge(std::make_pair("1", "3"));
  diGraph.addEdge(std::make_pair("3", "2"));
  diGraph.addEdge(std::make_pair("3", "4"));
  ASSERT_EQ(diGraph.getSize(), 4);

  std::string v1("1"), v2("3");
  auto const subgraphVertices = std::vector<gsl::not_null<std::string const *>>({&v1, &v2});
  auto const subgraph = graph.getSubgraph(subgraphVertices);
  ASSERT_EQ(subgraph->getSize(), 1);
  ASSERT_TRUE(subgraph->hasVertex(v1));
  ASSERT_TRUE(subgraph->hasVertex(v2));

  ASSERT_TRUE(subgraph->hasEdge(std::make_pair(&v1, &v2)));
  ASSERT_TRUE(subgraph->hasEdge(std::make_pair(&v2, &v1)));

  auto const subDigraph = diGraph.getSubgraph(subgraphVertices);
  ASSERT_EQ(subgraph->getSize(), 1);
  ASSERT_TRUE(subDigraph->hasVertex(v1));
  ASSERT_TRUE(subDigraph->hasVertex(v2));

  ASSERT_TRUE(subDigraph->hasEdge(std::make_pair(&v1, &v2)));
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

  graph.addEdge(std::make_pair("1", "2"));
  graph.addEdge(std::make_pair("2", "3"));
  graph.addEdge(std::make_pair("2", "4"));
  graph.addEdge(std::make_pair("3", "4"));
  graph.addEdge(std::make_pair("4", "1"));
  ASSERT_EQ(graph.getSize(), 5);

  auto diGraph = lazybastard::graph::DiGraph();
  diGraph.addVertex(std::move(spVertex1));
  diGraph.addVertex(std::move(spVertex2));
  diGraph.addVertex(std::move(spVertex3));
  diGraph.addVertex(std::move(spVertex4));
  ASSERT_EQ(diGraph.getOrder(), 4);

  diGraph.addEdge(std::make_pair("1", "2"));
  diGraph.addEdge(std::make_pair("2", "3"));
  diGraph.addEdge(std::make_pair("2", "4"));
  diGraph.addEdge(std::make_pair("3", "4"));
  diGraph.addEdge(std::make_pair("4", "1"));
  ASSERT_EQ(diGraph.getSize(), 5);

  std::string v1("1"), v2("4");
  auto shortestPathVertices = std::make_pair(lazybastard::util::make_not_null_and_const(graph.getVertex(v1).get()),
                                             lazybastard::util::make_not_null_and_const(graph.getVertex(v2).get()));
  auto shortestPath =
      lazybastard::GraphUtil::getShortestPath(lazybastard::util::make_not_null_and_const(&graph), shortestPathVertices);
  ASSERT_EQ(shortestPath.size(), 2);
  ASSERT_EQ(shortestPath[0], "1");
  ASSERT_EQ(shortestPath[1], "4");

  shortestPathVertices = std::make_pair(lazybastard::util::make_not_null_and_const(diGraph.getVertex(v1).get()),
                                        lazybastard::util::make_not_null_and_const(diGraph.getVertex(v2).get()));
  shortestPath = lazybastard::GraphUtil::getShortestPath(lazybastard::util::make_not_null_and_const(&diGraph),
                                                         shortestPathVertices);
  ASSERT_EQ(shortestPath.size(), 3);
  ASSERT_EQ(shortestPath[0], "1");
  ASSERT_EQ(shortestPath[1], "2");
  ASSERT_EQ(shortestPath[2], "4");
}
#include <gsl/pointers>
#include <memory>

#include <gtest/gtest.h>

#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

TEST(GraphTest, BasicTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::make_pair(spVertex1.get(), spVertex2.get())), "1,2");

  auto vertexIDs = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto graph = lazybastard::graph::Graph();
  graph.addVertex(std::move(spVertex1));
  graph.addVertex(std::move(spVertex2));
  ASSERT_EQ(graph.getOrder(), 2);

  auto edgeID = graph.addEdge(vertexIDs);
  ASSERT_EQ(edgeID, "1,2");

  auto const v1 = std::string("1");
  auto const v2 = std::string("2");
  auto expected = std::make_pair(gsl::make_not_null(&v1), gsl::make_not_null(&v2));
  ASSERT_TRUE(graph.hasEdge(expected));
}

TEST(GraphTest, EdgeCreationTest) {
  auto const spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto const spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);

  auto const edge = lazybastard::graph::Edge(std::make_pair(spVertex1, spVertex2));

  auto const vertices = edge.getVertices();
  ASSERT_EQ(vertices.first->getID(), "1");
  ASSERT_EQ(vertices.second->getID(), "2");
}

TEST(GraphTest, InverseOrderTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  ASSERT_EQ(lazybastard::graph::Edge::getEdgeID(std::make_pair(spVertex1.get(), spVertex2.get())), "1,2");

  auto const vertexIDs = std::make_pair(spVertex1->getID(), spVertex2->getID());
  auto graph = lazybastard::graph::Graph();
  graph.addVertex(std::move(spVertex1));
  graph.addVertex(std::move(spVertex2));
  ASSERT_EQ(graph.getOrder(), 2);

  auto const edgeID = graph.addEdge(vertexIDs);
  ASSERT_EQ(edgeID, "1,2");

  auto const v1 = std::string("1");
  auto const v2 = std::string("2");
  auto expected = std::make_pair(gsl::make_not_null(&v1), gsl::make_not_null(&v2));
  ASSERT_TRUE(graph.hasEdge(expected));
}
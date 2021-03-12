#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "Prokrastinator.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

TEST(CCTest, BasicTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>("1", 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>("2", 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>("3", 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>("4", 0);
  auto spVertex5 = std::make_shared<lazybastard::graph::Vertex>("5", 0);

  auto graph = lazybastard::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  graph.addVertex(spVertex4->getSharedPtr());
  graph.addVertex(spVertex5->getSharedPtr());

  graph.addEdge(std::make_pair(spVertex1->getID(), spVertex2->getID()));
  graph.addEdge(std::make_pair(spVertex1->getID(), spVertex5->getID()));
  graph.addEdge(std::make_pair(spVertex2->getID(), spVertex3->getID()));
  graph.addEdge(std::make_pair(spVertex2->getID(), spVertex5->getID()));
  graph.addEdge(std::make_pair(spVertex3->getID(), spVertex5->getID()));
  graph.addEdge(std::make_pair(spVertex3->getID(), spVertex4->getID()));
  graph.addEdge(std::make_pair(spVertex4->getID(), spVertex5->getID()));

  auto cc = lazybastard::getConnectedComponents(&graph);
  ASSERT_EQ(cc.size(), 1);

  auto const contains = [](std::vector<lazybastard::graph::Vertex *> &v, lazybastard::graph::Vertex *pS) {
    return std::find(std::begin(v), std::end(v), pS) != std::end(v);
  };

  ASSERT_TRUE(contains(cc[0], spVertex1.get()));
  ASSERT_TRUE(contains(cc[0], spVertex2.get()));
  ASSERT_TRUE(contains(cc[0], spVertex3.get()));
  ASSERT_TRUE(contains(cc[0], spVertex4.get()));
  //  ASSERT_TRUE(contains(cc[0], spVertex5.get()));

  auto spVertex6 = std::make_shared<lazybastard::graph::Vertex>("6", 0);
  auto spVertex7 = std::make_shared<lazybastard::graph::Vertex>("7", 0);
  auto spVertex8 = std::make_shared<lazybastard::graph::Vertex>("8", 0);

  graph.addVertex(spVertex6->getSharedPtr());
  graph.addVertex(spVertex7->getSharedPtr());
  graph.addVertex(spVertex8->getSharedPtr());

  graph.addEdge(std::make_pair(spVertex6->getID(), spVertex7->getID()));
  graph.addEdge(std::make_pair(spVertex6->getID(), spVertex8->getID()));
  graph.addEdge(std::make_pair(spVertex7->getID(), spVertex8->getID()));

  cc = lazybastard::getConnectedComponents(&graph);
  ASSERT_EQ(cc.size(), 2);

  std::size_t firstIdx = 0, secondIdx = 1;
  if (cc[0].size() < cc[1].size()) {
    std::swap(firstIdx, secondIdx);
  }

  ASSERT_TRUE(contains(cc[firstIdx], spVertex1.get()));
  ASSERT_TRUE(contains(cc[firstIdx], spVertex2.get()));
  ASSERT_TRUE(contains(cc[firstIdx], spVertex3.get()));
  ASSERT_TRUE(contains(cc[firstIdx], spVertex4.get()));
  ASSERT_TRUE(contains(cc[firstIdx], spVertex5.get()));
  ASSERT_TRUE(contains(cc[secondIdx], spVertex6.get()));
  ASSERT_TRUE(contains(cc[secondIdx], spVertex7.get()));
  ASSERT_TRUE(contains(cc[secondIdx], spVertex8.get()));
}
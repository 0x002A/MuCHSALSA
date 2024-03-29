#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "Kernel.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"
#include "types/Direction.h"

TEST(CCTest, BasicTest) {
  auto spVertex1 = std::make_shared<muchsalsa::graph::Vertex>(1, 0);
  auto spVertex2 = std::make_shared<muchsalsa::graph::Vertex>(2, 0);
  auto spVertex3 = std::make_shared<muchsalsa::graph::Vertex>(3, 0);
  auto spVertex4 = std::make_shared<muchsalsa::graph::Vertex>(4, 0);
  auto spVertex5 = std::make_shared<muchsalsa::graph::Vertex>(5, 0);

  auto graph = muchsalsa::graph::Graph();
  graph.addVertex(spVertex1->getSharedPtr());
  graph.addVertex(spVertex2->getSharedPtr());
  graph.addVertex(spVertex3->getSharedPtr());
  graph.addVertex(spVertex4->getSharedPtr());
  graph.addVertex(spVertex5->getSharedPtr());

  graph.addEdge(std::make_pair(spVertex1.get(), spVertex2.get()));
  graph.addEdge(std::make_pair(spVertex1.get(), spVertex5.get()));
  graph.addEdge(std::make_pair(spVertex2.get(), spVertex3.get()));
  graph.addEdge(std::make_pair(spVertex2.get(), spVertex5.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex5.get()));
  graph.addEdge(std::make_pair(spVertex3.get(), spVertex4.get()));
  graph.addEdge(std::make_pair(spVertex4.get(), spVertex5.get()));

  // Set Consensus Direction
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex2.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex5.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex3.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex5.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex5.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex4.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex4.get(), spVertex5.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);

  auto cc = muchsalsa::getConnectedComponents(graph);
  ASSERT_EQ(cc.size(), 1);

  auto const contains = [](std::vector<muchsalsa::graph::Vertex *> &v, muchsalsa::graph::Vertex *pS) {
    return std::find(std::begin(v), std::end(v), pS) != std::end(v);
  };

  ASSERT_TRUE(contains(cc[0], spVertex1.get()));
  ASSERT_TRUE(contains(cc[0], spVertex2.get()));
  ASSERT_TRUE(contains(cc[0], spVertex3.get()));
  ASSERT_TRUE(contains(cc[0], spVertex4.get()));
  ASSERT_TRUE(contains(cc[0], spVertex5.get()));

  auto spVertex6 = std::make_shared<muchsalsa::graph::Vertex>(6, 0);
  auto spVertex7 = std::make_shared<muchsalsa::graph::Vertex>(7, 0);
  auto spVertex8 = std::make_shared<muchsalsa::graph::Vertex>(8, 0);

  graph.addVertex(spVertex6->getSharedPtr());
  graph.addVertex(spVertex7->getSharedPtr());
  graph.addVertex(spVertex8->getSharedPtr());

  graph.addEdge(std::make_pair(spVertex6.get(), spVertex7.get()));
  graph.addEdge(std::make_pair(spVertex6.get(), spVertex8.get()));
  graph.addEdge(std::make_pair(spVertex7.get(), spVertex8.get()));

  graph.getEdge(std::make_pair(spVertex6.get(), spVertex7.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex6.get(), spVertex8.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);
  graph.getEdge(std::make_pair(spVertex7.get(), spVertex8.get()))->setConsensusDirection(muchsalsa::Direction::e_POS);

  cc = muchsalsa::getConnectedComponents(graph);
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

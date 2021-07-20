#include <gtest/gtest.h>

#include "Prokrastinator.h"
#include "graph/Edge.h"
#include "graph/Graph.h"
#include "graph/Vertex.h"

TEST(MSTTest, BasicTest) {
  auto spVertex1 = std::make_shared<lazybastard::graph::Vertex>(1, 0);
  auto spVertex2 = std::make_shared<lazybastard::graph::Vertex>(2, 0);
  auto spVertex3 = std::make_shared<lazybastard::graph::Vertex>(3, 0);
  auto spVertex4 = std::make_shared<lazybastard::graph::Vertex>(4, 0);
  auto spVertex5 = std::make_shared<lazybastard::graph::Vertex>(5, 0);

  auto graph = lazybastard::graph::Graph();
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

  // Set weights
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex2.get()))->setWeight(4);
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex5.get()))->setWeight(7);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex3.get()))->setWeight(3);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex5.get()))->setWeight(5);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex5.get()))->setWeight(6);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex4.get()))->setWeight(3);
  graph.getEdge(std::make_pair(spVertex4.get(), spVertex5.get()))->setWeight(1);

  auto mst = lazybastard::getMaxSpanTree(graph);
  ASSERT_EQ(mst.getSize(), 0);

  // Set consensus direction
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex2.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex1.get(), spVertex5.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex3.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex2.get(), spVertex5.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex5.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex3.get(), spVertex4.get()))->setConsensusDirection(true);
  graph.getEdge(std::make_pair(spVertex4.get(), spVertex5.get()))->setConsensusDirection(true);

  mst = lazybastard::getMaxSpanTree(graph);
  ASSERT_EQ(mst.getSize(), 4);
  ASSERT_EQ(mst.getOrder(), graph.getOrder());

  ASSERT_TRUE(mst.hasEdge(std::make_pair(spVertex1.get(), spVertex5.get())));
  ASSERT_TRUE(mst.hasEdge(std::make_pair(spVertex3.get(), spVertex4.get())));
  ASSERT_TRUE(mst.hasEdge(std::make_pair(spVertex5.get(), spVertex2.get())));
  ASSERT_TRUE(mst.hasEdge(std::make_pair(spVertex5.get(), spVertex3.get())));
  ASSERT_FALSE(mst.hasEdge(std::make_pair(spVertex1.get(), spVertex2.get())));
  ASSERT_FALSE(mst.hasEdge(std::make_pair(spVertex2.get(), spVertex3.get())));
  ASSERT_FALSE(mst.hasEdge(std::make_pair(spVertex4.get(), spVertex5.get())));
}

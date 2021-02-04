#pragma once

namespace lazybastard {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace coroutine {

template <typename T> class generator;

} // namespace coroutine

namespace graph {

class DiGraph;
class Edge;
class EdgeOrder;
class Graph;
class Vertex;

} // namespace graph

namespace matching {

class MatchMap;
class ID2OverlapMap;

} // namespace matching

namespace threading {

class Job;
class ThreadPool;

} // namespace threading

namespace util {

template <typename T> struct LTCmp;

} // namespace util

class OutputWriter;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

} // namespace lazybastard
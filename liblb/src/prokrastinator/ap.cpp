#include "Prokrastinator.h"

#include "OutputWriter.h"
#include "graph/Graph.h"
#include "matching/ID2OverlapMap.h"
#include "matching/MatchMap.h"

void lazybastard::assemblePath(gsl::not_null<graph::Graph const *> const /*pGraph*/,
                               gsl::not_null<lazybastard::matching::MatchMap const *> const /*matches*/,
                               gsl::not_null<lazybastard::matching::ID2OverlapMap *> const /*pID2OverlapMap*/,
                               gsl::not_null<std::vector<gsl::not_null<std::string const *>> const *> const /*pPath*/,
                               gsl::not_null<lazybastard::graph::DiGraph const *> const /*pDiGraph*/,
                               std::size_t /*idx*/, lazybastard::OutputWriter & /*writer*/) {}
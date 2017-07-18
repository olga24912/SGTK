#include "RuleDelSmallCycle.h"

namespace filter {
    namespace scaffolder {
        void filter::scaffolder::RuleDelSmallCycle::simplifyGraph(filter::ContigGraph *filter) {
            std::vector<int> vert = filter->getVertexList();
            for (int v : vert) {
                std::vector<int> edges = filter->getEdges(v);
                std::vector<int> egdesR = filter->getEdgesR(v);

                for (int e : edges) {
                    for (int er : egdesR) {
                        int u = filter->getEdgeTo(e);
                        int w = filter->getEdgeFrom(er);
                        if (u == w) {
                            if (MAX_DIF * filter->getEdgeWeight(e) <= filter->getEdgeWeight(er)) {
                                std::stringstream ss;
                                ss << e;
                                filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                            } else if (MAX_DIF * filter->getEdgeWeight(er) <= filter->getEdgeWeight(e)) {
                                std::stringstream ss;
                                ss << er;
                                filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                            } else {
                                std::stringstream ss;
                                ss << e;
                                filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                                ss.clear();
                                ss << er;
                                filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                            }
                        }
                    }
                }
            }
        }
    }
}

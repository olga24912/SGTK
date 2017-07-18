#include "RuleBigDifInWeight.h"

namespace filter {
    namespace scaffolder {
        void RuleBigDifInWeight::simplifyGraph(filter::ContigGraph *filter) {
            INFO("start simplify graph");
            std::vector<int> vert = filter->getVertexList();
            for (int v : vert) {
                delSmallEdges(filter, filter->getEdges(v));
                delSmallEdges(filter, filter->getEdgesR(v));
            }
            INFO("finish simplify graph");
        }

        void RuleBigDifInWeight::delSmallEdges(ContigGraph *filter, const std::vector<int> &edges) const {
            int maxW = 0;
            for (int e : edges) {
                if (filter->getEdgeWeight(e) >= maxW) {
                    maxW = filter->getEdgeWeight(e);
                }
            }

            for (int e : edges) {
                if (filter->getEdgeWeight(e) * maxDif <= maxW) {
                    std::stringstream ss;
                    ss << e;
                    filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                }
            }
        }
    }
}

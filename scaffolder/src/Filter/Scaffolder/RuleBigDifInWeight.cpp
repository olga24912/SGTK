#include "RuleBigDifInWeight.h"

namespace filter {
    namespace scaffolder {
        void RuleBigDifInWeight::simplifyGraph(filter::Filter *filter) {
            std::vector<int> vert = filter->getVertexList();
            for (int v : vert) {
                delSmallEdges(filter, filter->getEdges(v));
                delSmallEdges(filter, filter->getEdgesR(v));
            }
        }

        void RuleBigDifInWeight::delSmallEdges(Filter *filter, const std::vector<int> &edges) const {
            int maxW = 0;
            for (int e : edges) {
                if (filter->getEdgeWeight(e) >= maxW) {
                    maxW = filter->getEdgeWeight(e);
                }
            }

            for (int e : edges) {
                if (filter->getEdgeWeight(e) * maxDif < maxW) {
                    std::stringstream ss;
                    ss << e;
                    filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                }
            }
        }
    }
}

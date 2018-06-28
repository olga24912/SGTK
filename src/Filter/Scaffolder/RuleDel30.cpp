#include "RuleDel30.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;

        void RuleDel30::simplifyGraph(ContigGraph *filter) {
            std::vector<int> vert = filter->getVertexList();

            for (int v: vert) {
                std::vector<int> edges = filter->getEdges(v);

                for (int e : edges) {
                    if (filter->getLibType(filter->getEdgeLib(e)) == ContigGraph::Lib::RNA_SPLIT_30) {
                        filter->delEdge(e);
                    }
                }
            }
        }
    }
}

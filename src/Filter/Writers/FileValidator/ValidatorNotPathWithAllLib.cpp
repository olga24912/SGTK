#include <iostream>
#include "ValidatorNotPathWithAllLib.h"

namespace filter {
    namespace writers {
        bool ValidatorNotPathWithAllLib::isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) {
            DEBUG("start validation");
            std::vector<int> listLib = graph->getLibList();

            for (int i = 0; i < (int) vert.size(); ++i) {
                int v = vert[i];
                std::vector<int> edges = graph->getEdges(v);
                if (edges.size() == 0) continue;
                int u = graph->getEdgeTo(edges[0]);

                std::vector<int> was(listLib[listLib.size() - 1] + 1, 0);
                for (int e : edges) {
                    if (graph->getEdgeTo(e) != u) {
                        return true;
                    }
                    was[graph->getEdgeLib(e)] = 1;
                }

                for (int l : listLib) {
                    if (was[l] == 0) {
                        return true;
                    }
                }
            }

            return false;
        }
    }
}
#include <iostream>
#include "ValidatorNotPathWithAllLib.h"

bool ValidatorNotPathWithAllLib::isGoodVertexSet(std::vector<int> vert, Filter *filter) {
    std::cerr <<"start val" << std::endl;
    std::vector<int> listLib = filter->getLibList();

    for (int i = 0; i < (int)vert.size(); ++i) {
        int v = vert[i];
        std::vector<int> edges = filter->getEdges(v);
        if (edges.size() == 0) continue;
        int u = filter->getEdgeTo(edges[0]);

        std::vector<int> was(listLib[listLib.size() - 1], 0);
        for (int e : edges) {
            if (filter->getEdgeTo(e) != u) {
                return true;
            }
            was[filter->getEdgeLib(e)] = 1;
        }

        for (int l : listLib) {
            if (was[l] == 0) {
                return true;
            }
        }
    }

    return false;
}

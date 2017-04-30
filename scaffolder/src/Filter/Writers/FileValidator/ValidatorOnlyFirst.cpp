#include <algorithm>
#include "ValidatorOnlyFirst.h"

bool ValidatorOnlyFirst::isGoodVertexSet(std::vector<int> vert, Filter *filter) {
    for (int i = 0; i < (int)vert.size(); ++i) {
        int v = vert[i];

        std::vector<int> edges = filter->getEdges(v);
        std::vector<std::pair <int, int> > ed;

        for (int e : edges ) {
            int u = filter->getEdgeTo(e);
            int inSet = 0;
            for (int y : vert) {
                if (u == y) inSet = 1;
            }
            if (inSet) {
                ed.push_back(std::make_pair(filter->getEdgeTo(e), filter->getEdgeLib(e)));
            }
        }
        std::sort(ed.begin(), ed.end());
        ed.push_back(std::make_pair(-1, 0));

        int wasPr = 0, wasNotPr = 0;

        for (int j = 0; j < (int)ed.size(); ++j) {
            if (j > 0 && ed[j].first != ed[j - 1].first) {
                if (wasPr && !wasNotPr) {
                    return true;
                }
            }
            if (ed[j].second == libPr) wasPr = 1;
            if (ed[j].second == libNotPr) wasNotPr = 1;
        }
    }

    return false;
}

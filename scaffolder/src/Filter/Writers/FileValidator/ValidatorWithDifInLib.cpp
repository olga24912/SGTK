#include <iostream>
#include "ValidatorWithDifInLib.h"
#include <algorithm>

bool ValidatorWithDifInLib::isGoodVertexSet(std::vector<int> vert, Filter *filter) {
    std::cerr <<"start val" << std::endl;
    std::vector<int> listLib = filter->getLibList();


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

        std::vector<int> was((unsigned)listLib[listLib.size() - 1] + 1, 0);

        ed.push_back(std::make_pair(-1, 0));

        for (int j = 0; j < (int)ed.size(); ++j) {
            if (j > 0 && ed[j].first != ed[j - 1].first) {
                int cnt = 0;
                for (int l : libs) {
                    if (was[l] != 0) {
                        ++cnt;
                    }
                    was[l] = 0;
                }
                if (cnt > 0 && cnt < libs.size()) {
                    return true;
                }
            }
            was[ed[j].second] = 1;
        }
    }

    return false;
}

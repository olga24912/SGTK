#include <iostream>
#include "ValidatorWithDifInLib.h"
#include <algorithm>

namespace filter {
    namespace writers {
        bool ValidatorWithDifInLib::isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) {
            DEBUG("start validation");
            std::vector<int> listLib = graph->getLibList();


            for (int i = 0; i < (int) vert.size(); ++i) {
                int v = vert[i];

                std::vector<int> edges = graph->getEdges(v);
                std::vector<std::pair<int, int> > ed;

                for (int e : edges) {
                    int u = graph->getEdgeTo(e);
                    int inSet = 0;
                    for (int y : vert) {
                        if (u == y) inSet = 1;
                    }
                    if (inSet) {
                        ed.push_back(std::make_pair(graph->getEdgeTo(e), graph->getEdgeLib(e)));
                    }
                }
                std::sort(ed.begin(), ed.end());

                std::vector<int> was((unsigned) listLib[listLib.size() - 1] + 1, 0);

                ed.push_back(std::make_pair(-1, 0));

                for (int j = 0; j < (int) ed.size(); ++j) {
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
    }
}
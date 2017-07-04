#include <algorithm>
#include <set>
#include "ScaffoldStrategy.h"

namespace filter {
    namespace scaffolder {
        void ScaffoldStrategy::topSort(Filter *graph) {
            INFO("top sort");
            int n = graph->getVertexCount();
            topSortPos.resize(n);

            std::vector<int> used(n, 0);

            for (int i = 0; i < n; ++i) {
                if (used[i] == 0) {
                    topSortDfs(i, graph, &used);
                }
            }

            std::reverse(topsort.begin(), topsort.end());
            for (int i = 0; i < (int) topsort.size(); ++i) {
                topSortPos[topsort[i]] = i;
            }
        }

        void ScaffoldStrategy::findCycle(Filter *graph) {
            INFO("condensation");
            int n = graph->getVertexCount();
            color.resize(n);
            int col = 1;
            for (int i = (int) topsort.size() - 1; i >= 0; --i) {
                if (color[topsort[i]] == 0) {
                    colorDfs(topsort[i], col, graph);
                    ++col;
                }
            }
        }


        void ScaffoldStrategy::topSortDfs(int v, Filter *graph, std::vector<int> *used) {
            TRACE("topSortDfs v=" << v);
            (*used)[v] = 1;
            std::vector<int> edges = graph->getEdges(v);
            for (int e : edges) {
                int u = graph->getEdgeTo(e);
                if ((*used)[u] == 0) {
                    topSortDfs(u, graph, used);
                }
            }

            topsort.push_back(v);
        }

        void ScaffoldStrategy::colorDfs(int v, int col, Filter *graph) {
            TRACE("colorDfs v=" << v);
            color[v] = col;

            std::vector<int> edges = graph->getEdges(v);
            for (int e : edges) {
                int u = graph->getEdgeTo(e);
                if (color[u] == 0) {
                    colorDfs(u, col, graph);
                }
            }
        }


        int ScaffoldStrategy::deg(int i, Filter *pFilter, int dirIn) {
            std::vector<int> edges;
            if (dirIn) {
                edges = pFilter->getEdgesR(i);
            } else {
                edges = pFilter->getEdges(i);
            }

            std::set<int> nb;
            for (int e : edges) {
                if (dirIn) {
                    nb.insert(pFilter->getEdgeFrom(e));
                } else {
                    nb.insert(pFilter->getEdgeTo(e));
                }
            }

            TRACE("deg i=" << i << ": " << nb.size());
            return (int) nb.size();
        }
    }
}
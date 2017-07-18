#include "RuleDelCycle.h"

namespace filter {
    namespace scaffolder {
        void RuleDelCycle::simplifyGraph(filter::ContigGraph *filter) {
            INFO("start del cycles");
            int wasOpt = true;

            while(wasOpt) {
                wasOpt = false;
                topSort(filter);
                int cntCol = findCycle(filter);

                TRACE("Start find min edge in cycle =" << cntCol);
                std::vector<int> minColor(cntCol, 1);

                std::vector<int> vect = filter->getVertexList();
                /*for (int v : vect) {
                    std::vector<int> edges = graph->getEdges(v);

                    for (int e : edges) {
                        int u = graph->getEdgeTo(e);
                        if (color[v] == color[u]) {
                            if (minColor[color[v]] == 0 ||
                                    minColor[color[v]] > graph->getEdgeWeight(e)) {
                                minColor[color[v]] = graph->getEdgeWeight(e);
                            }
                        }
                    }
                }*/

                for (int v : vect) {
                    std::vector<int> edges = filter->getEdges(v);

                    for (int e : edges) {
                        int u = filter->getEdgeTo(e);
                        if (color[v] == color[u]) {
                            if (minColor[color[v]] * DEL_DIF >= filter->getEdgeWeight(e)) {
                                std::stringstream ss;
                                ss << e;
                                filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                                //wasOpt = true;
                            }
                        }
                    }
                }
            }
            INFO("Finish del cycles");
        }

        void  RuleDelCycle::topSort(ContigGraph *graph) {
            DEBUG("top sort");
            int n = graph->getVertexCount();
            topSortPos.resize(0);
            topSortPos.resize(n);
            topsort.resize(0);

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
            DEBUG("finish topsort");
        }

        int  RuleDelCycle::findCycle(ContigGraph *graph) {
            DEBUG("condensation");
            int n = graph->getVertexCount();
            color.resize(0);
            color.resize(n);
            int col = 1;
            for (int i = (int) topsort.size() - 1; i >= 0; --i) {
                if (color[topsort[i]] == 0) {
                    colorDfs(topsort[i], col, graph);
                    ++col;
                }
            }
            DEBUG("finish condensation");
            return col;
        }


        void  RuleDelCycle::topSortDfs(int v, ContigGraph *graph, std::vector<int> *used) {
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
            TRACE("finish topsort dfs v=" << v);
        }

        void  RuleDelCycle::colorDfs(int v, int col, ContigGraph *graph) {
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
    }
}

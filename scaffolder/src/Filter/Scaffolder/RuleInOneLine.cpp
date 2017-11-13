#include <queue>
#include <set>
#include "RuleInOneLine.h"

namespace filter {
    namespace scaffolder {
        void RuleInOneLine::simplifyGraph(ContigGraph *graph) {
            std::vector<int> vect = graph->getVertexList();
            for (int v : vect) {
                std::vector<int> edges = graph->getEdges(v);

                if (edges.size() > 20) continue; //TODO skip big deg

                for (int i = 0; i < (int) edges.size(); ++i) {
                    projectEdge(graph, edges[i]);
                }
            }
        }

        void RuleInOneLine::projectEdge(ContigGraph *graph, int e) {
            int v = graph->getEdgeFrom(e);
            int w = graph->getEdgeTo(e);

            int cnt_path = 0;

            std::vector<int> path;
            std::vector<int> edges = reduceEdges(graph, graph->getEdges(v), graph->getLibType(graph->getEdgeLib(e)));
            for (int e2 : edges) {
                if (graph->getEdgeTo(e2) != w) {
                    int u = graph->getEdgeTo(e2);

                    int cnt = getPaths(graph, u, w, path, graph->getLibType(graph->getEdgeLib(e)));
                    cnt_path += cnt;
                    if (cnt == 1) {
                        path.push_back(e2);
                    }
                    if (cnt_path > 1) return;
                }
            }
            if (cnt_path == 0) return;

            for (int e2 : path) {
                if (graph->getLibType(graph->getEdgeLib(e2)) == graph->getLibType(graph->getEdgeLib(e))) {
                    graph->setWeight(e2, graph->getEdgeWeight(e2) + graph->getEdgeWeight(e));
                } else {
                    graph->addEdge(graph->getEdgeFrom(e2), graph->getEdgeTo(e2),
                                   graph->getEdgeLib(e), graph->getEdgeWeight(e), graph->getEdgeCoordB1(e2),
                    graph->getEdgeCoordE1(e2), graph->getEdgeCoordB2(e2), graph->getEdgeCoordE2(e2));
                }
            }
            graph->delEdge(e);
        }

        int RuleInOneLine::getPaths(ContigGraph *graph, int u, int w, std::vector<int> &path, ContigGraph::Lib::Type type) {
            DEBUG("start get path");
            std::set<int> uarea = getArea(graph, u);
            std::set<int> warea = getArea(graph, w);

            DEBUG("uarea size =" << uarea.size() << " warea size = " <<warea.size() );

            if (warea.count(u) > 0 && uarea.count(w) > 0) return 2;
            if (uarea.count(w) == 0) return 0;

            std::vector<int> topSort = topSortF(graph, uarea);
            std::set<int> cycle = markCycle(topSort, graph, uarea);

            std::map<int, int> d;
            d[u] = 1;
            std::map<int, int> p;
            p[u] = -1;

            //DEBUG("topsort: ")
            for (int i = 0; i < (int) topSort.size(); ++i) {
                //DEBUG(topSort[i]);
                int v = topSort[i];
                std::vector<int> edges = reduceEdges(graph, graph->getEdges(v), type);

                for (int e : edges) {
                    d[graph->getEdgeTo(e)] += d[v];
                    p[graph->getEdgeTo(e)] = e;
                }
            }

            if (d[w] != 1) return d[w];

            DEBUG("start find path");

            int cur = w;
            while (p[cur] != -1) {
                if (cycle.count(cur)) return 2;
                path.push_back(p[cur]);
                cur = graph->getEdgeFrom(p[cur]);
                //DEBUG("cur = " << cur << " u=" << u << " w=" << w << " d[cur] = " << d[cur] << " d[u] =" << d[u] << " d[w]=" << d[w] << " inCycle:" << cycle.count(cur));
            }

            return 1;
        }

        std::set<int> RuleInOneLine::getArea(ContigGraph *graph, int v) {
            std::unordered_map<int, int> dist;
            std::set<int> elems;
            std::queue<int> que;
            que.push(v);
            dist[v] = 0;
            elems.insert(v);

            while (que.size() != 0) {
                int u = que.front();
                que.pop();

                if (dist[v] == DIST || elems.size() > 50) continue;

                std::vector<int> edges = graph->getEdges(u);
                for (int edge : edges) {
                    int y = graph->getEdgeTo(edge);
                    if (dist.count(y) == 0) {
                        dist[y] = dist[u] + 1;
                        que.push(y);
                        elems.insert(y);
                    }
                }
            }

            return elems;
        }

        std::vector<int> RuleInOneLine::topSortF(ContigGraph *graph, std::set<int> uarea) {
            std::set<int> used;
            std::vector<int> topSort;
            for (int v : uarea) {
                if (used.count(v) == 0) {
                    topSortDfs(v, graph, uarea, used, topSort);
                }
            }

            std::reverse(topSort.begin(), topSort.end());

            return topSort;
        }

        void RuleInOneLine::topSortDfs(int v, ContigGraph *graph, std::set<int> &area,
                                       std::set<int> &used, std::vector<int> &topSort) {
            used.insert(v);
            std::vector<int> edges = graph->getEdges(v);
            for (int e : edges) {
                int u = graph->getEdgeTo(e);
                if (area.count(u) == 0) continue;
                if (used.count(u) == 1) continue;

                topSortDfs(u, graph, area, used, topSort);
            }

            topSort.push_back(v);
        }

        std::set<int> RuleInOneLine::markCycle(std::vector<int> &topSort, ContigGraph *graph, std::set<int> area) {
            std::set<int> used;
            std::set<int> cycle;
            for (int i = topSort.size() - 1; i >= 0; --i) {
                if (used.count(topSort[i]) == 0) {
                    dfsMarkCycle(topSort[i], graph, area, used, cycle, false);
                }
            }

            return cycle;
        }

        void RuleInOneLine::dfsMarkCycle(int &v, ContigGraph *graph, std::set<int> area, std::set<int>& used,
                                         std::set<int>& cycle, bool inCycle) {
            if (inCycle) {
                cycle.insert(v);
            }

            std::vector<int> edges = graph->getEdges(v);
            used.insert(v);

            for (int e : edges) {
                int u = graph->getEdgeTo(e);
                if (area.count(u) == 0) continue;
                if (used.count(u) == 1) continue;

                cycle.insert(v);
                dfsMarkCycle(u, graph, area, used, cycle, true);
            }
        }

        std::vector<int>
        RuleInOneLine::reduceEdges(ContigGraph *graph, std::vector<int> edges, ContigGraph::Lib::Type type) {
            std::vector<int> res;
            for (int e : edges) {
                int was = 0;
                for (int i = 0; i < (int)res.size(); ++i) {
                    int e1 = res[i];
                    if (sameEdges(graph, e, e1)) {
                        if (graph->getLibType(graph->getEdgeLib(e)) == type) {
                            res[i] = e;
                        }
                        was = 1;
                    }
                }
                if (was == 0) {
                    res.push_back(e);
                }
            }

            return res;
        }

        bool RuleInOneLine::sameEdges(ContigGraph *graph, int e1, int e2) {
            return (graph->getEdgeFrom(e1) == graph->getEdgeFrom(e2)) &&
                    (graph->getEdgeTo(e1) == graph->getEdgeTo(e2)) &&
                    (sameCoord1(graph, e1, e2) && sameCoord2(graph, e1, e2));
        }
    }
}

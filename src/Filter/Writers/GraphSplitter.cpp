#include <queue>
#include <iostream>
#include "GraphSplitter.h"

namespace filter {
    namespace writers {
        std::vector<std::vector<int> > GraphSplitter::split(ContigGraph *graph1, std::vector<int> vert) {
            INFO("start split file on parts");
            for (int v : vert) {
                std::cerr << v << " ";
            }
            std::cerr << "\n";

            DEBUG("vert size = " << vert.size());

            this->graph = graph1;
            this->vert = vert;

            clear();

            for (int v : vert) {
                while (used[v] == 0) {
                    findNewComp(v);
                }
            }

            INFO("finish split file on parts");
            return res;
        }

        void GraphSplitter::clear() {
            DEBUG("clear");
            used.resize(0);
            res.resize(0);
            edgeCol.resize(0);

            used.resize((unsigned int) graph->getVertexCount(), 2);
            for (int v : vert) {
                if (v >= used.size()) {
                    used.resize(v + 1, 2);
                }
                used[v] = 0;
            }


            int edgeCnt = 0;
            for (int i = 0; i < (unsigned) used.size(); ++i) {
                for (int e : graph->getEdges(i)) {
                    edgeCnt = std::max(edgeCnt, e + 1);
                }
            }

            edgeCol.resize(edgeCnt);
            DEBUG("finish clear");
        }

        void GraphSplitter::findNewComp(int v) {
            TRACE("find new components for v=" << v << " edges=" << graph->getEdges(v).size());
            std::queue<int> que;
            int cntV = 0;
            int cntE = 0;

            int colNum = (int) res.size();
            res.push_back(std::vector<int>());

            que.push(v);
            while (que.size() > 0) {
                if (cntV == maxVert) break;
                int u = que.front();
                que.pop();

                int extraEdge = 0;
                for (int e : graph->getEdges(u)) {
                    int w = graph->getEdgeTo(e);
                    if (used[w] == 1) {
                        ++extraEdge;
                    }
                }

                for (int e : graph->getEdgesR(u)) {
                    int w = graph->getEdgeFrom(e);
                    if (used[w] == 1) {
                        ++extraEdge;
                    }
                }

                if (cntE + extraEdge <= maxEdge && cntV + 1 <= maxVert) {
                    res[colNum].push_back(u);
                    cntE += extraEdge;
                    cntV += 1;
                    used[u] = 1;
                    for (int e : graph->getEdges(u)) {
                        int w = graph->getEdgeTo(e);
                        if (used[w] == 1) {
                            edgeCol[e] = colNum + 1;
                        } else if (edgeCol[e] == 0 && used[w] == 0) {
                            used[w] = 3;
                            que.push(w);
                        }
                    }

                    for (int e : graph->getEdgesR(u)) {
                        int w = graph->getEdgeFrom(e);
                        if (used[w] == 1) {
                            edgeCol[e] = colNum + 1;
                        } else if (edgeCol[e] == 0 && used[w] == 0) {
                            used[w] = 3;
                            que.push(w);
                        }
                    }

                } else {
                    used[u] = 3;
                }
            }

            for (int i = 0; i < (int) used.size(); ++i) {
                if (used[i] == 3) used[i] = 0;
                if (used[i] == 1) {
                    used[i] = 2;
                    for (int e : graph->getEdges(i)) {
                        if (edgeCol[e] == 0 && used[graph->getEdgeTo(e)] != 2) {
                            used[i] = 0;
                        }
                    }


                    for (int e : graph->getEdgesR(i)) {
                        if (edgeCol[e] == 0 && used[graph->getEdgeFrom(e)] != 2) {
                            used[i] = 0;
                        }
                    }
                }
            }
        }
    }
}
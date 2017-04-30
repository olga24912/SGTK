#include <algorithm>
#include <iostream>
#include <set>
#include "ScaffolderPipeline.h"

void ScaffolderPipeline::evaluate(Filter *graph, std::string out) {

    topSort(graph);
    findCycle(graph);
    ignoreLightEdge(graph);
    uniqueConnection(graph);
    inOneLineConnection(graph);
    addConnectionForLeaves(graph, 0);
    addConnectionForLeaves(graph, 1);

    //oneWightConnection(graph);

    scaffolds.print(out);
}

void ScaffolderPipeline::uniqueConnection(Filter* graph) {
    for (int i = 0; i < (int)topsort.size(); ++i) {
        int v = topsort[i];
        std::vector<int> edges = (graph->getEdges(v));
        if (edges.size() == 0) {
            continue;
        }
        int u = graph->getEdgeTo(edges[0]);
        if (isUniquePair(v, u, graph) && color[v] != color[u]) {
            scaffolds.addConnection(v, u);
        }
    }
}

bool ScaffolderPipeline::isUniquePair(int v1, int v2, Filter * graph) {
    std::vector<int> edges = graph->getEdges(v1);
    for (int e : edges) {
        int u = graph->getEdgeTo(e);
        if (u != v2) {
            return false;
        }
    }

    edges = graph->getEdgesR(v2);
    for (int e : edges) {
        int u = graph->getEdgeFrom(e);
        if (u != v1) {
            return false;
        }
    }

    return true;
}

void ScaffolderPipeline::topSort(Filter *graph) {
    int n = graph->getVertexCount();
    topSortPos.resize(n);

    std::vector<int> used(n, 0);

    for (int i = 0; i < n; ++i) {
        if (used[i] == 0) {
            topSortDfs(i, graph, &used);
        }
    }

    std::reverse(topsort.begin(), topsort.end());
    for (int i = 0; i < (int)topsort.size(); ++i) {
        topSortPos[topsort[i]] = i;
    }
}

void ScaffolderPipeline::findCycle(Filter *graph) {
    int n = graph->getVertexCount();
    color.resize(n);
    int col = 1;
    for (int i = (int)topsort.size() - 1; i >= 0; --i) {
        if (color[topsort[i]] == 0) {
            colorDfs(topsort[i], col, graph);
            ++col;
        }
    }
}

void ScaffolderPipeline::topSortDfs(int v, Filter *graph, std::vector<int>* used) {
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

void ScaffolderPipeline::colorDfs(int v, int col, Filter * graph) {
    color[v] = col;

    std::vector<int> edges = graph->getEdges(v);
    for (int e : edges) {
        int u = graph->getEdgeTo(e);
        if (color[u] == 0) {
            colorDfs(u, col, graph);
        }
    }
}

ScaffolderPipeline::ScaffolderPipeline(std::string contigFile) : scaffolds(Scaffolds(contigFile)) {
}

void ScaffolderPipeline::inOneLineConnection(Filter *graph) {
    addFirstConnection(graph);
    delEdgeFromDifPath(graph);
}

void ScaffolderPipeline::addFirstConnection(Filter *graph) {
    for (int i = 0; i < (int)topsort.size(); ++i) {
        int v = topsort[i];

        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;
        int minu = graph->getEdgeTo(edges[0]);

        for (int e : edges) {
            int u = graph->getEdgeTo(e);

            if (topSortPos[u] < topSortPos[minu] && color[v] != color[u]) {
                minu = u;
            }
        }

        if (color[v] != color[minu]) {
            scaffolds.addConnection(v, minu);
        }
    }
}

void ScaffolderPipeline::delEdgeFromDifPath(Filter *graph) {
    for (int i = (int)topsort.size() - 1; i >= 0; --i) {
        int v = topsort[i];
        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;

        for (int e : edges) {
            int u = graph->getEdgeTo(e);
            assert(graph->getEdgeFrom(e) == v);

            if (scaffolds.lineId(v) != scaffolds.lineId(u)) {
                scaffolds.brokeConnection(v);
                scaffolds.brokeConnectionTo(u);
            }

        }
    }
}

void ScaffolderPipeline::oneWightConnection(Filter *graph) {
    std::vector<int> libs = graph->getLibList();
    for (int l : libs) {
        if (l == 0 || l >= 6) continue;
        std::stringstream ss;
        ss << l << " " << 1;
        graph->processQuery(Query(Query::MIN_EDGE_WEIGHT, ss.str()));
    }

    uniqueConnection(graph);
}

void ScaffolderPipeline::addConnectionForLeaves(Filter *graph, int dir) {
    int n = graph->getVertexCount();

    for (int i = 0; i < n; ++i) {
        std::vector<int> edges = dir == 0 ? graph->getEdges(i) : graph->getEdgesR(i);


        int u = -1, v = -1;
        for (int e : edges) {
            int to = dir == 0 ? graph->getEdgeTo(e) : graph->getEdgeFrom(e);
            if (u == -1 || u == to) {
                u = to;
            } else if (v == -1 || v == to) {
                v = to;
            } else {
                u = -2;
                break;
            }
        }
        if (i == 3187) {
            std::cerr <<"3187 "<<  u << " " << v << std::endl;
        }

        if (u < 0 || v < 0) {
            continue;
        }

        if (color[i] == color[u] || color[i] == color[v] || color[u] == color[v]) {
            continue;
        }

        if (deg(v, graph, dir) != 0) {
            std::swap(v, u);
        }
        if (deg(i, graph, !dir) <= 1 &&
                deg(v, graph, !dir) == 1 && deg(u, graph, !dir) == 1 &&
                deg(v, graph, dir) == 0 && deg(u, graph, dir) <= 1) {
            if (dir == 0) {
                scaffolds.addConnection(i, u);
            } else {
                scaffolds.addConnection(u, i);
            }
        }

    }
}

int ScaffolderPipeline::deg(int i, Filter *pFilter, int dirIn) {
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

    return (int) nb.size();
}

void ScaffolderPipeline::ignoreLightEdge(Filter *graph) {
    std::vector<int> vert = graph->getVertexList();
    for (int v : vert) {
        std::vector<int> maxWInLib(graph->getLibList()[graph->getLibList().size() - 1], 0);

        std::vector<int> edges = graph->getEdges(v);
        for (int e : edges) {
            maxWInLib[graph->getEdgeLib(e)] =
                    std::max(maxWInLib[graph->getEdgeLib(e)], graph->getEdgeWeight(e));
        }

        for (int e : edges) {
            if (graph->getEdgeWeight(e) * 2 < maxWInLib[graph->getEdgeLib(e)]) {
                std::stringstream ss;
                ss << e;
                graph->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
            }
        }
    }
}
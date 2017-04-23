#include <algorithm>
#include <iostream>
#include "ScaffolderPipeline.h"

void ScaffolderPipeline::evaluate(Filter *graph, std::string out) {

    topSort(graph);
    findCycle(graph);
    uniqueConnection(graph);

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
        if (v == 659) {
            std::cerr << v << " " << u << " " << isUniquePair(v, u, graph)<< " " << color[v] << " " << color[u] << std::endl;
        }
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

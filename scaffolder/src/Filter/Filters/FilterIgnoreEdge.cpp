#include "FilterIgnoreEdge.h"

FilterIgnoreEdge::FilterIgnoreEdge(Filter *filter) : Filter(filter) {
}

void FilterIgnoreEdge::processQuery(Query query) {
    if (query.type == query.SET_IGNORE_EDGE) {
        std::stringstream ss(query.argv);
        int e;
        ss >> e;
        ignore.resize(std::max(e + 1, ignore.size()), 0);
        ignore[e] = 1;
        e ^= 1;
        ignore.resize(std::max(e + 1, ignore.size()), 0);
        ignore[e] = 1;
    } else {
        Filter::processQuery(query);
    }
}

std::vector<int> FilterIgnoreEdge::getEdges(int v) {
    std::vector<int> edges = subfilter->getEdges(v);
    std::vector<int> res;

    for (int e : edges) {
        if (e >= ignore.size() || ignore[e] == 0) {
            res.push_back(e);
        }
    }

    return res;
}

std::vector<int> FilterIgnoreEdge::getEdgesR(int v) {
    std::vector<int> edges = subfilter->getEdgesR(v);
    std::vector<int> res;

    for (int e : edges) {
        if (e >= ignore.size() || ignore[e] == 0) {
            res.push_back(e);
        }
    }

    return res;
}

#include "FilterIgnore.h"

std::vector<int> FilterIgnore::getVertexList() {
    std::vector<int> vertList = subfilter->getVertexList();
    std::vector<int> res;
    for (int i = 0; i < (int)vertList.size(); ++i) {
        int v = vertList[i];
        if (ignore[v] == 0) {
            res.push_back(v);
        }
    }
    return res;
}

std::vector<int> FilterIgnore::getEdges(int v) {
    std::vector<int> edges = subfilter->getEdges(v);
    std::vector<int> res;
    for (int e : edges) {
        int u = subfilter->getEdgeTo(e);
        if (ignore[u] == 0) {
            res.push_back(e);
        }
    }

    return res;
}

std::vector<int> FilterIgnore::getEdgesR(int v) {
    std::vector<int> edges = subfilter->getEdgesR(v);
    std::vector<int> res;
    for (int e : edges) {
        int u = subfilter->getEdgeFrom(e);
        if (ignore[u] == 0) {
            res.push_back(e);
        }
    }

    return res;
}

FilterIgnore::FilterIgnore(Filter *filter) : Filter(filter) {
    subfilter = filter;
    int cnt = subfilter->getVertexCount();

    ignore.resize(cnt, 0);
}

void FilterIgnore::processQuery(Query query) {
    if (query.type == query.SET_IGNORE) {
        std::stringstream ss(query.argv);
        int vstart, vend;
        ss >> vstart >> vend;
        for (int i = vstart; i < vend; ++i) {
            ignore[i] = 1;
        }
    } else if (query.type == query.RESET_IGNORE) {
        for (int i = 0; i < (int)ignore.size(); ++i) {
            ignore[i] = 0;
        }
    } else {
        Filter::processQuery(query);

        ignore.resize(subfilter->getVertexCount(), 0);
    }
}

#include <iostream>
#include "FilterMinWeight.h"

FilterMinWeight::FilterMinWeight(Filter *filter) : Filter(filter) {
    if ((int)filter->getLibList().size() > 0) {
        libMinEdgeWight.resize((unsigned) filter->getLibList()
        [(int) filter->getLibList().size() - 1] + 1, 1);
    } else {
        libMinEdgeWight.resize(0);
    }
}

std::vector<int> FilterMinWeight::getEdges(int v) {
    std::vector<int> edges = subfilter->getEdges(v);
    std::vector<int> res;
    for (int i = 0; i < (int)edges.size(); ++i) {
        int e = edges[i];
        if (libMinEdgeWight[subfilter->getEdgeLib(e)] <= subfilter->getEdgeWeight(e)) {
            int u = subfilter->getEdgeTo(e);
            if (minContigLen <= subfilter->getTargetLen(u)) {
                res.push_back(e);
            }
        }
    }

    return res;
}

std::vector<int> FilterMinWeight::getEdgesR(int v) {
    std::vector<int> edges = subfilter->getEdgesR(v);
    std::vector<int> res;
    for (int i = 0; i < (int)edges.size(); ++i) {
        int e = edges[i];
        if (libMinEdgeWight[subfilter->getEdgeLib(e)] <= subfilter->getEdgeWeight(e)) {
            int v = subfilter->getEdgeFrom(e);
            if (minContigLen <= subfilter->getTargetLen(v)) {
                res.push_back(e);
            }
        }
    }

    return res;
}

std::vector<int> FilterMinWeight::getVertexList() {
    std::vector<int> vertList = subfilter->getVertexList();
    std::vector<int> res;
    for (int i = 0; i < (int)vertList.size(); ++i) {
        int v = vertList[i];
        if (minContigLen <= subfilter->getTargetLen(v)) {
            res.push_back(v);
        }
    }

    return res;
}

void FilterMinWeight::processQuery(Query query) {
    if (query.type == query.MIN_EDGE_WEIGHT) {
        std::stringstream ss(query.argv);
        int libNum;
        int weight;
        ss >> libNum >> weight;
        libMinEdgeWight[libNum] = weight;
    } else if (query.type == query.MIN_CONTIG_LEN){
        std::stringstream ss(query.argv);
        int weight;
        ss >> weight;
        minContigLen = weight;
    } else {
        Filter::processQuery(query);
        if (subfilter->getLibList().size() > 0) {
            libMinEdgeWight.resize((unsigned) subfilter->getLibList()
            [subfilter->getLibList().size() - 1] + 1, 0);
        }
    }
}

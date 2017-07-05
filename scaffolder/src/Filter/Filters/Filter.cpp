#include <fstream>
#include <iostream>
#include <algorithm>
#include "Filter.h"

namespace filter {
    int Filter::getVertexCount() {
        assert(subfilter != nullptr);
        TRACE("getVertexCount");
        return subfilter->getVertexCount();
    }

    std::vector<int> Filter::getEdges(int v) {
        assert(subfilter != nullptr);
        TRACE("getEdges v=" << v);
        return subfilter->getEdges(v);
    }

    std::vector<int> Filter::getEdgesR(int v) {
        assert(subfilter != nullptr);
        TRACE("getEdgesR v=" << v);
        return subfilter->getEdgesR(v);
    }

    std::string Filter::getTargetName(int v) {
        assert(subfilter != nullptr);
        TRACE("getTargetName v=" << v);
        return subfilter->getTargetName(v);
    }

    int Filter::getTargetLen(int v) {
        assert(subfilter != nullptr);
        TRACE("getTargetLen v=" << v);
        return subfilter->getTargetLen(v);
    }

    int Filter::getEdgeTo(int e) {
        assert(subfilter != nullptr);
        TRACE("getEdgeTo e=" << e);
        return subfilter->getEdgeTo(e);
    }

    int Filter::getEdgeFrom(int e) {
        assert(subfilter != nullptr);
        TRACE("getEdgeFrom e=" << e);
        return subfilter->getEdgeFrom(e);
    }

    int Filter::getEdgeWeight(int e) {
        assert(subfilter != nullptr);
        TRACE("getEdgeWeight e=" << e);
        return subfilter->getEdgeWeight(e);
    }

    std::vector<int> Filter::getLibList() {
        assert(subfilter != nullptr);
        TRACE("getLibList");
        return subfilter->getLibList();
    }

    int Filter::getEdgeLib(int e) {
        assert(subfilter != nullptr);
        TRACE("getEdgeLib e=" << e);
        return subfilter->getEdgeLib(e);
    }

    std::string Filter::getLibName(int l) {
        assert(subfilter != nullptr);
        TRACE("getLibName l=" << l);
        return subfilter->getLibName(l);
    }

    std::string Filter::getLibColor(int l) {
        assert(subfilter != nullptr);
        TRACE("getLibColor l=" << l);
        return subfilter->getLibColor(l);
    }

    std::vector<int> Filter::getVertexList() {
        assert(subfilter != nullptr);
        TRACE("getVertexList");
        return subfilter->getVertexList();
    }

    void Filter::processQuery(Query query) {
        assert(subfilter != nullptr);
        TRACE("processQuery");
        return subfilter->processQuery(query);
    }

    void Filter::write(std::string fileName) {
        INFO("write graph to fileName=" << fileName);
        std::ofstream out(fileName);

        std::vector<int> libList = getLibList();
        out << libList.size() << "\n";
        std::vector<int> newColr(libList[libList.size() - 1] + 1, 0);
        for (int i = 0; i < (int) libList.size(); ++i) {
            out << "l " << i << " " << getLibColor(libList[i]) << " " << getLibName(libList[i]) << "\n";
            newColr[libList[i]] = i;
        }

        std::vector<int> edges;
        out << getVertexCount() << "\n";
        for (int i = 0; i < getVertexCount(); ++i) {
            out << "v " << i << " " << getTargetName(i) << " " << getTargetLen(i) << "\n";
            std::vector<int> curEd = getEdges(i);
            for (int j = 0; j < (int) curEd.size(); ++j) {
                edges.push_back(curEd[j]);
            }
        }

        std::sort(edges.begin(), edges.end());
        edges.resize(unique(edges.begin(), edges.end()) - edges.begin());
        out << edges.size() << "\n";
        for (int i = 0; i < edges.size(); ++i) {
            out << "e " << i << " " << getEdgeFrom(edges[i]) << " " << getEdgeTo(edges[i]) << " "
                << newColr[getEdgeLib(edges[i])] << " " << getEdgeWeight(edges[i]) << "\n";
        }

        out.close();
    }

    std::string Filter::getInfo(int e) {
        TRACE("getInfo e=" << e);
        return subfilter->getInfo(e);
    }


    contig_graph::ContigGraph::Lib::Type Filter::getLibType(int l) {
        TRACE("get lib type l=" << l);
        return subfilter->getLibType(l);
    }

    int Filter::getEdgeCoordB1(int e) {
        TRACE("getEdgeCoord begin1 e=" << e);
        return subfilter->getEdgeCoordB1(e);
    }

    int Filter::getEdgeCoordE1(int e) {
        TRACE("getEdgeCoord end1 e=" << e);
        return subfilter->getEdgeCoordE1(e);
    }

    int Filter::getEdgeCoordB2(int e) {
        TRACE("get Edge coords begin2 e=" << e);
        return subfilter->getEdgeCoordB2(e);
    }

    int Filter::getEdgeCoordE2(int e) {
        TRACE("get edge coords ebd2 e=" << e);
        return subfilter->getEdgeCoordE2(e);
    }
}
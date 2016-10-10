//
// Created by olga on 08.10.16.
//

#include "DNAPairReadGraphBuilder.h"

void DNAPairReadGraphBuilder::setDistBetweenPairReads(int distBetweenPairReads) {
    DNAPairReadGraphBuilder::distBetweenPairReads = distBetweenPairReads;
}

void DNAPairReadGraphBuilder::addInfoAboutRead(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAboutRead(readName, target, read);
    read1DistToEnd[readName] = readDist(read);
}

int DNAPairReadGraphBuilder::readDist(BamAlignmentRecord read) {
    if (!hasFlagRC(read)) {
        return (graph->getTargetLength(2 * read.rID) - read.beginPos);
    } else {
        return (read.beginPos + read.tLen);
    }
}

void DNAPairReadGraphBuilder::addInfoAbout2Read(string readName, int target, BamAlignmentRecord read) {
    PairReadGraphBuilder::addInfoAbout2Read(readName, target, read);
    read2DistToEnd[readName] = readDist(read);
}

void DNAPairReadGraphBuilder::filterEdge() {
    PairReadGraphBuilder::filterEdge();

    for (int v = 0; v < graph->getVertexCount(); ++v) {
        graph->sortEdgeByWeight(v);
        vector<int> edges_weight = graph->getEdgesWeight(v);
        if (edges_weight.size() == 0) {
            continue;
        }
        sort(edges_weight.begin(), edges_weight.end());
        int cnt = countEdgesBeforeBreak(v, edges_weight);
        graph->delEdges(v, (int)edges_weight.size() - cnt);
    }
}

int DNAPairReadGraphBuilder::countEdgesBeforeBreak(int v, vector<int> edges) {
    int cnt = 1;
    int max_val = 0;
    if (edges.size() > 0) {
        max_val = edges[edges.size() - 1];
    }
    for (int i = static_cast<int>(edges.size()) - 2; i >= 0; --i) {
        int w0 = edges[i + 1], w1 = edges[i];
        if ((w0 - w1) < max_val * DEFAULT_DEF) {
            ++cnt;
        } else {
            break;
        }
    }
    if (cnt > DEFAULT_MAX_COUNT_EDGE) {
        cnt = 0;
    }

    return cnt;
}

void DNAPairReadGraphBuilder::incEdgeWeight(string readName, int target) {
    if (read1Target.count(readName)) {
        if (read1Target[readName] == target ||
            read1Target[readName] == pairTarget(target)) {
            return;
        }

        if (read1DistToEnd[readName] + read2DistToEnd[readName] >
            distBetweenPairReads) {
            return;
        }

        int verFID = read1Target[readName], verSID = target,
                verRFID = pairTarget(verFID), verRSID = pairTarget(verSID);
        graph->incEdgeWeight(verFID, verSID);
        graph->incEdgeWeight(verRSID, verRFID);
    }
}

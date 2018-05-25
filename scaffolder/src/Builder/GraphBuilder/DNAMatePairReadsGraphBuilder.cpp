#include <tgmath.h>
#include "DNAMatePairReadsGraphBuilder.h"


int builder::graph_builder::DNAMatePairReadsGraphBuilder::get1Target(const seqan::BamAlignmentRecord &read) const {
    TRACE("get1Target");

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev != firstRev) {
        target++;
    }
    return target;
}


int builder::graph_builder::DNAMatePairReadsGraphBuilder::get2Target(const seqan::BamAlignmentRecord &read) const {
    TRACE("get2Target");

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev != secondRev) {
        target++;
    }
    return target;
}

void builder::graph_builder::DNAMatePairReadsGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord& read1,
                                                                 seqan::BamAlignmentRecord& read2) {
    TRACE("incEdgeWeight read1 " << read1.beginPos << " "
                                 << (read1.beginPos + seqan::getAlignmentLengthInRef(read1)) << " RC=" <<  hasFlagRC(read1) <<
                                 " target " << get1Target(read1));
    TRACE("incEdgeWeight read2 " << read2.beginPos << " "
                                 << (read2.beginPos + seqan::getAlignmentLengthInRef(read2)) << " RC=" << hasFlagRC(read2) <<
                                 " target " << get2Target(read2));

    int target1 = get1Target(read1);
    std::pair<int, int> t1c = std::make_pair(0, 0);
    int target2 = get2Target(read2);
    std::pair<int, int> t2c = std::make_pair(0, 0);

    if (target1 < 0 || hasFlagSecondary(read1) || !isUniqueMapRead(read1)) {
        return;
    }
    if (target2 < 0 || hasFlagSecondary(read2) || !isUniqueMapRead(read2)) {
        return;
    }

    if (target1 == target2 || target1 == pairTarget(target2)) {
        return;
    }

    changeEdges(target1, t1c, target2, t2c);
    changeEdges(pairTarget(target2),
                std::make_pair(graph->getTargetLen(target2) - t2c.second, graph->getTargetLen(target2) - t2c.first),
                pairTarget(target1),
                std::make_pair(graph->getTargetLen(target1) - t1c.second, graph->getTargetLen(target1) - t1c.first));

}

int builder::graph_builder::DNAMatePairReadsGraphBuilder::changeEdges(int v1, std::pair<int, int> c1, int v2,
                                                              std::pair<int, int> c2) {
    std::vector<ContigGraph::Edge> edges = graph->getEdgesBetween(v1, v2);

    for (ContigGraph::Edge edge: edges) {
        graph->incEdge(edge.id, c1, c2);
        return edge.id;
    }

    int e = graph->addEdge(v1, v2, c1, c2);
    return e;
}
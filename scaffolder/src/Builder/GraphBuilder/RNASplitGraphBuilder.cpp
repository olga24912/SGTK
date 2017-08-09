#include "RNASplitGraphBuilder.h"

int builder::graph_builder::RNASplitGraphBuilder::get1Target(const seqan::BamAlignmentRecord &read) const {
    TRACE("get1Target");

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev) {
        target++;
    }
    return target;
}


int builder::graph_builder::RNASplitGraphBuilder::get2Target(const seqan::BamAlignmentRecord &read) const {
    TRACE("get2Target");

    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev) {
        target++;
    }
    return target;
}

void builder::graph_builder::RNASplitGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord read1,
                                                                 seqan::BamAlignmentRecord read2) {
    INFO("incEdgeWeight read1 " << read1.beginPos << " "
                                << (read1.beginPos + seqan::getAlignmentLengthInRef(read1)) << " RC=" <<  hasFlagRC(read1) <<
    " target " << get1Target(read1));
    INFO("incEdgeWeight read2 " << read2.beginPos << " "
                                << (read2.beginPos + seqan::getAlignmentLengthInRef(read2)) << " RC=" << hasFlagRC(read2) <<
    " target " << get2Target(read2));

    assert(seqan::isUniqueMapRead(read1));
    assert(seqan::isUniqueMapRead(read2));

    int target1 = get1Target(read1);
    std::pair<int, int> t1c = getCoord(read1, target1);
    int target2 = get2Target(read2);
    std::pair<int, int> t2c = getCoord(read2, target2);

    if (target1 == target2 || target1 == pairTarget(target2)) {
        return;
    }

    int e1 = changeEdges(target1, t1c, target2, t2c);
    int e2 = changeEdges(pairTarget(target2),
                std::make_pair(graph->getTargetLen(target2) - t2c.second, graph->getTargetLen(target2) - t2c.first),
                pairTarget(target1),
                std::make_pair(graph->getTargetLen(target1) - t1c.second, graph->getTargetLen(target1) - t1c.first));

    samFileWriter.writeEdge(e1, read1, read2);
    samFileWriter.writeEdge(e2, read2, read1);
}

std::pair<int, int> builder::graph_builder::RNASplitGraphBuilder::getCoord(seqan::BamAlignmentRecord read, int target) {
    if ((hasFlagRC(read))) {
        return std::make_pair(graph->getTargetLen(target) - (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)),
                              graph->getTargetLen(target) - read.beginPos);
    } else {
        return std::make_pair(read.beginPos, (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)));
    }
}

int builder::graph_builder::RNASplitGraphBuilder::changeEdges(int v1, std::pair<int, int> c1, int v2,
                                                              std::pair<int, int> c2) {
    std::vector<ContigGraph::Edge> edges = graph->getEdgesBetween(v1, v2);

    for (ContigGraph::Edge edge: edges) {
        if (isGoodEdgeFor1(edge, c1) && isGoodEdgeFor2(edge, c2)) {
            std::pair<int, int> coord1 = relaxCoord(std::make_pair(edge.coordBegin1, edge.coordEnd1), c1);
            std::pair<int, int> coord2 = relaxCoord(std::make_pair(edge.coordBegin2, edge.coordEnd2), c2);

            graph->incEdge(edge.id, coord1, coord2);

            return edge.id;
        }
    }

    int e = graph->addEdge(v1, v2, c1, c2);
    return e;
}

bool builder::graph_builder::RNASplitGraphBuilder::isGoodEdgeFor1(builder::contig_graph::ContigGraph::Edge edge, std::pair<int, int> c) {
    if (c.first < c.second) {
        if (edge.coordBegin1 >= edge.coordEnd1) return false;

        return std::abs(edge.coordEnd1 - c.second) < 10;
    } else {
        if (edge.coordBegin1 <= edge.coordEnd1) return false;

        return std::abs(edge.coordBegin1 - c.first) < 10;
    }
}

bool builder::graph_builder::RNASplitGraphBuilder::isGoodEdgeFor2(builder::contig_graph::ContigGraph::Edge edge,
                                                                  std::pair<int, int> c) {
    if (c.first < c.second) {
        if (edge.coordBegin1 >= edge.coordEnd1) return false;

        return std::abs(edge.coordBegin1 - c.first) < 10;
    } else {
        if (edge.coordBegin1 <= edge.coordEnd1) return false;

        return std::abs(edge.coordEnd1 - c.second) < 10;
    }
}

std::pair<int, int>
builder::graph_builder::RNASplitGraphBuilder::relaxCoord(std::pair<int, int> c1,
                                                         std::pair<int, int> c) {
    if (c.first < c.second) {
        return std::make_pair(std::min(c.first, c1.first), std::max(c.second, c1.second));
    } else {
        return std::make_pair(std::max(c.first, c1.first), std::min(c.second, c1.second));
    }
}

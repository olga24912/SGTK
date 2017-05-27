#include <cmath>
#include "BlockSplitDotWriter.h"

void BlockSplitDotWriter::writeOneVertexSet(std::vector<int> vert, std::string fileName) {
    std::vector<std::vector<vertBlock> > blocks(vert.size());
    for (int i = 0; i < (int)vert.size(); ++i) {
        blocks[i] = splitOnBlocks(vert[i]);
    }

    findOutsideEdge(blocks);

    for (int i = 0; i < (int)blocks.size(); ++i) {
        writeOneVertBlock(blocks[i]);
    }

    for (int i = 0; i < (int)blocks.size(); ++i) {
        writeBlockEdges(blocks[i]);
    }

    writeEdges(blocks);
}

BlockSplitDotWriter::BlockSplitDotWriter(Filter *filter, FileValidator *validator, int maxVert, int maxEdge)
        : DotWriter(filter, validator, maxVert, maxEdge) {}

std::vector<BlockSplitDotWriter::vertBlock> BlockSplitDotWriter::splitOnBlocks(int v) {
    std::vector<vertBlock> blocks;
    std::vector<int> edges = filter->getEdges(v);

    for (int e : edges) {
        int cb = filter->getEdgeCoordB1(e), ce = filter->getEdgeCoordE1(e);

        int was = 0;

        for (vertBlock b : blocks) {
            if (std::fabs(b.coordE - cb) < minBlockDist || std::fabs(b.coordB - ce) < minBlockDist) {
                b.coordB = std::min(b.coordB, cb);
                b.coordE = std::min(b.coordE, ce);
                was = 1;
            }
        }

        blocks.push_back(vertBlock(v, cb, ce));
    }

    edges = filter->getEdgesR(v);

    for (int e : edges) {
        int cb = filter->getEdgeCoordB2(e), ce = filter->getEdgeCoordE2(e);

        for (vertBlock b : blocks) {
            if (std::fabs(b.coordE - cb) < minBlockDist || std::fabs(b.coordB - ce) < minBlockDist) {
                b.coordB = std::min(b.coordB, cb);
                b.coordE = std::min(b.coordE, ce);
            }
        }
        blocks.push_back(vertBlock(v, cb, ce));
    }

    if (filter->getTargetName(v)[filter->getTargetName(v).size() - 1] == 'v') {
        std::sort(blocks.begin(), blocks.end());
    } else {
        std::sort(blocks.rbegin(), blocks.rend());
    }

    return blocks;
}

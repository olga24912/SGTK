#include <cmath>
#include "BlockSplitDotWriter.h"

void BlockSplitDotWriter::writeOneVertexSet(std::vector<int> vert, std::string fileName) {
    std::vector<std::vector<vertBlock> > blocks(vert.size());
    for (int i = 0; i < (int)vert.size(); ++i) {
        blocks[i] = splitOnBlocks(vert[i]);
    }


    findOutsideEdge(blocks);

    std::ofstream out(fileName);
    out << "digraph {\n";
    for (int i = 0; i < (int)blocks.size(); ++i) {
        writeOneVertBlock(blocks[i], out);
    }
    for (int i = 0; i < (int)blocks.size(); ++i) {
        writeBlockEdges(blocks[i]);
    }
    writeEdges(blocks);

    out << "}\n";
    out.close();
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

        if (was == 0) blocks.push_back(vertBlock(v, cb, ce));
    }

    edges = filter->getEdgesR(v);

    for (int e : edges) {
        int cb = filter->getEdgeCoordB2(e), ce = filter->getEdgeCoordE2(e);
        int was = 0;

        for (vertBlock b : blocks) {
            if (std::fabs(b.coordE - cb) < minBlockDist || std::fabs(b.coordB - ce) < minBlockDist) {
                b.coordB = std::min(b.coordB, cb);
                b.coordE = std::min(b.coordE, ce);
                was = 1;
            }
        }
        if (was == 0) blocks.push_back(vertBlock(v, cb, ce));
    }

    if (filter->getTargetName(v)[filter->getTargetName(v).size() - 1] == 'v') {
        std::sort(blocks.begin(), blocks.end());
    } else {
        std::sort(blocks.rbegin(), blocks.rend());
    }

    return blocks;
}

void BlockSplitDotWriter::findOutsideEdge(const std::vector<std::vector<BlockSplitDotWriter::vertBlock>>& blocks) {
    for (int i = 0; i < (int)blocks.size(); ++i) {
        std::vector<int> edges = filter->getEdges(blocks[i][0].vertId);

        for (int e : edges) {
            if (isOutsideEdge(e, blocks)) {
                for (int j = 0; j < (int)blocks[i].size(); ++j) {
                    if (blocks[i][j].coordB <= filter->getEdgeCoordB1(e) &&
                            blocks[i][j].coordE >= filter->getEdgeCoordE1(e)) {
                        blocks[i][j].hasOutsideEdge = 1;
                    }
                }
            }
        }


        edges = filter->getEdgesR(blocks[i][0].vertId);

        for (int e : edges) {
            if (isOutsideEdge(e, blocks)) {
                for (int j = 0; j < (int)blocks[i].size(); ++j) {
                    if (blocks[i][j].coordB <= filter->getEdgeCoordB2(e) &&
                        blocks[i][j].coordE >= filter->getEdgeCoordE2(e)) {
                        blocks[i][j].hasOutsideEdge = 1;
                    }
                }
            }
        }

    }
}

bool BlockSplitDotWriter::isOutsideEdge(int e, const std::vector<std::vector<BlockSplitDotWriter::vertBlock>> &bl) {
    int v = filter->getEdgeFrom(e);
    int u = filter->getEdgeTo(e);

    int was1 = 0, was2 = 0;

    for (int i = 0; i < (int)bl.size(); ++i) {
        if (bl[i][0].vertId == v) {
            was1 = 1;
        }
        if (bl[i][0].vertId == u) {
            was2 = 1;
        }
    }

    return was1 == 0 || was2 == 0;
}

void BlockSplitDotWriter::writeOneVertBlock(std::vector<BlockSplitDotWriter::vertBlock> &bl, std::ofstream &out) {
    for (int i = 0; i < (int)bl.size(); ++i) {
        std::stringstream ss;
        ss << filter->getTargetName(bl[i].vertId) << "_" << i;
        int v = bl[i].vertId;

        out << "    \"" << ss.str() << "\"[label=\" " << filter->getTargetName(v) <<
            " id = " << v
            << "\nlen = " << filter->getTargetLen(v) <<
            << "\ncoord:" << bl[i].coordB << "-"<< bl[i].coordE << "\"";
        if (bl[i].hasOutsideEdge) {
            out << " , style = \"filled\", color = \"#F0E68C\"";
        }
        out << "];\n";
    }
}

void BlockSplitDotWriter::writeBlockEdges(const std::vector<BlockSplitDotWriter::vertBlock> &bl, std::ofstream &out) {
    for (int i = 0; i < (int)bl.size() - 1; ++i) {
        std::stringstream ssname1;
        std::stringstream ssname2;
        int v = bl[i].vertId;
        ssname1 << filter->getTargetName(v) << "_" << i;
        ssname2 << filter->getTargetName(v) << "_" << i + 1;

        out << "    \"" << ssname1.str() << "\" -> \"";
        out << ssname2.str() << "\" [ ";
        out << "penwidth = " << 100 << "]\n";
    }
}
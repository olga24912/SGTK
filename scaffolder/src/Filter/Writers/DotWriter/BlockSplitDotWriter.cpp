#include <cmath>
#include "BlockSplitDotWriter.h"

namespace filter {
    namespace writers {
        void BlockSplitDotWriter::writeOneVertexSet(std::vector<int> vert, std::string fileName) {
            DEBUG("writeOneVertrexSet");

            if (vert.size() <= 1) return;
            std::vector<std::vector<vertBlock> > blocks(vert.size());
            for (int i = 0; i < (int) vert.size(); ++i) {
                blocks[i] = splitOnBlocks(vert[i]);
            }


            findOutsideEdge(blocks);

            std::ofstream out(fileName);
            out << "digraph {\n";
            out << "node[shape=box]\n";
            for (int i = 0; i < (int) blocks.size(); ++i) {
                writeOneVertBlock(blocks[i], out);
            }
            for (int i = 0; i < (int) blocks.size(); ++i) {
                writeBlockEdges(blocks[i], out);
            }
            writeEdges(blocks, out);

            out << "}\n";
            out.close();
        }

        BlockSplitDotWriter::BlockSplitDotWriter(ContigGraph *graph1, FileValidator *validator, int maxVert, int maxEdge)
                : DotWriter(graph1, validator, maxVert, maxEdge) {}

        std::vector<BlockSplitDotWriter::vertBlock> BlockSplitDotWriter::splitOnBlocks(int v) {
            TRACE("split on blocks v=" << v);
            std::vector<vertBlock> blocks;
            std::vector<int> edges = graph->getEdges(v);

            for (int e : edges) {
                int cb = graph->getEdgeCoordB1(e), ce = graph->getEdgeCoordE1(e);

                int was = 0;

                for (vertBlock b : blocks) {
                    if (std::fabs(b.coordE - cb) < minBlockDist || std::fabs(b.coordB - ce) < minBlockDist ||
                        (b.coordB <= cb && cb <= b.coordE) || (b.coordB <= ce && ce <= b.coordE) ||
                        (cb <= b.coordB && b.coordB <= ce)) {
                        b.coordB = std::min(b.coordB, cb);
                        b.coordE = std::max(b.coordE, ce);
                        was = 1;
                    }
                }

                if (was == 0) blocks.push_back(vertBlock(v, cb, ce));
            }

            edges = graph->getEdgesR(v);

            for (int e : edges) {
                int cb = graph->getEdgeCoordB2(e), ce = graph->getEdgeCoordE2(e);
                int was = 0;

                for (vertBlock b : blocks) {
                    if (std::fabs(b.coordE - cb) < minBlockDist || std::fabs(b.coordB - ce) < minBlockDist ||
                        (b.coordB <= cb && cb <= b.coordE) || (b.coordB <= ce && ce <= b.coordE) ||
                        (cb <= b.coordB && b.coordB <= ce)) {
                        b.coordB = std::min(b.coordB, cb);
                        b.coordE = std::max(b.coordE, ce);
                        was = 1;
                    }
                }
                if (was == 0) blocks.push_back(vertBlock(v, cb, ce));
            }

            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                std::sort(blocks.begin(), blocks.end());
            } else {
                std::sort(blocks.rbegin(), blocks.rend());
            }

            return blocks;
        }

        void BlockSplitDotWriter::findOutsideEdge(std::vector<std::vector<BlockSplitDotWriter::vertBlock>> &blocks) {
            TRACE("find outside edge for blocks");
            for (int i = 0; i < (int) blocks.size(); ++i) {

                TRACE("block size =" << blocks[i].size());
                if (blocks[i].size() == 0) continue;
                std::vector<int> edges = graph->getEdges(blocks[i][0].vertId);

                for (int e : edges) {
                    if (isOutsideEdge(e, blocks)) {
                        for (int j = 0; j < (int) blocks[i].size(); ++j) {
                            if (blocks[i][j].coordB <= graph->getEdgeCoordB1(e) &&
                                blocks[i][j].coordE >= graph->getEdgeCoordE1(e)) {
                                blocks[i][j].hasOutsideEdge = 1;
                            }
                        }
                    }
                }


                edges = graph->getEdgesR(blocks[i][0].vertId);

                for (int e : edges) {
                    if (isOutsideEdge(e, blocks)) {
                        for (int j = 0; j < (int) blocks[i].size(); ++j) {
                            if (blocks[i][j].coordB <= graph->getEdgeCoordB2(e) &&
                                blocks[i][j].coordE >= graph->getEdgeCoordE2(e)) {
                                blocks[i][j].hasOutsideEdge = 1;
                            }
                        }
                    }
                }

            }
        }

        bool
        BlockSplitDotWriter::isOutsideEdge(int e, const std::vector<std::vector<BlockSplitDotWriter::vertBlock>> &bl) {
            TRACE("isOutSizeEdge e=" << e);
            int v = graph->getEdgeFrom(e);
            int u = graph->getEdgeTo(e);

            int was1 = 0, was2 = 0;

            for (int i = 0; i < (int) bl.size(); ++i) {
                if (bl[i][0].vertId == v) {
                    was1 = 1;
                }
                if (bl[i][0].vertId == u) {
                    was2 = 1;
                }
            }

            return was1 == 0 || was2 == 0;
        }

        void
        BlockSplitDotWriter::writeOneVertBlock(std::vector<BlockSplitDotWriter::vertBlock> &bl, std::ofstream &out) {
            TRACE("writeOneVertBlock");
            std::stringstream ss;
            ss << "node" << bl[0].vertId << "_" << 0;
            int v = bl[0].vertId;
            out << "    \"" << ss.str() << "\"[label=\" " << graph->getTargetName(v) <<
                " id = " << v
                << "\nlen = " << graph->getTargetLen(v)
                << "\ncoord:" << bl[0].coordB << "-" << bl[0].coordE << "\"";
            if (bl[0].hasOutsideEdge) {
                out << " , style = \"filled\", color = \"#F0E68C\"";
            }
            out << "];\n";

            for (int i = 1; i < (int) bl.size(); ++i) {
                std::stringstream ss;
                ss << "node" << bl[i].vertId << "_" << i;
                int v = bl[i].vertId;

                out << "    \"" << ss.str() << "\"[label=\" coord:" << bl[i].coordB << "-" << bl[i].coordE << "\"";

                if (bl[i].hasOutsideEdge) {
                    out << " , style = \"filled\", color = \"#F0E68C\"";
                }
                out << "];\n";
            }

            out << "{rank=same;";
            for (int i = 0; i < (int) bl.size(); ++i) {
                out << "node" << bl[i].vertId << "_" << i << " ";

            }
            out << "}\n";
        }

        void
        BlockSplitDotWriter::writeBlockEdges(const std::vector<BlockSplitDotWriter::vertBlock> &bl,
                                             std::ofstream &out) {
            TRACE("write block edges");
            for (int i = 0; i < (int) bl.size() - 1; ++i) {
                std::stringstream ssname1;
                std::stringstream ssname2;
                int v = bl[i].vertId;
                ssname1 << "node" << v << "_" << i;
                ssname2 << "node" << v << "_" << i + 1;

                out << "    " << ssname1.str() << " -> ";
                out << ssname2.str() << " [ ";
                out << "penwidth = " << 10 << "]\n";
            }
        }

        void BlockSplitDotWriter::writeEdges(const std::vector<std::vector<BlockSplitDotWriter::vertBlock>> &bl,
                                             std::ofstream &out) {
            DEBUG("write edges");

            for (int i = 0; i < (int) bl.size(); ++i) {
                if (bl[i].size() == 0) continue;
                int v = bl[i][0].vertId;

                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    int u = graph->getEdgeTo(e);

                    int i1 = findBlockId(v, e, bl);
                    int i2 = findBlockId(u, e, bl);

                    if (i1 == -1 || i2 == -1) continue;
                    if (v == u) continue;


                    TRACE("v=" << v << " u=" << u << " e=" << e << " weight=" << graph->getEdgeWeight(e));
                    assert(graph->getEdgeWeight(e) > 0);

                    out << "    " << "node" << v << "_" << i1 << " -> ";
                    out << "node" << u << "_" << i2 << " [ ";
                    out << "color = \"" << graph->getLibColor(graph->getEdgeLib(e)) << "\", ";
                    out << "penwidth = " << 1 + (int) log10((double) graph->getEdgeWeight(e)) << ", ";
                    out << "label = " << "\"" << graph->getLibName(graph->getEdgeLib(e));
                    out << "\n weight = " << (graph->getEdgeWeight(e));
                    out << "\n id = " << e << "\" ]\n";
                }
            }
        }

        int
        BlockSplitDotWriter::findBlockId(int v, int e,
                                         const std::vector<std::vector<BlockSplitDotWriter::vertBlock>> &bl) {
            TRACE("find block id for v=" << v << " e=" << e);
            for (int i = 0; i < (int) bl.size(); ++i) {
                if (bl[i][0].vertId == v) {
                    for (int j = 0; j < (int) bl[i].size(); ++j) {
                        if (graph->getEdgeFrom(e) == v) {
                            if (graph->getEdgeCoordB1(e) >= bl[i][j].coordB &&
                                graph->getEdgeCoordE1(e) <= bl[i][j].coordE) {
                                return j;
                            }
                        } else {
                            if (graph->getEdgeCoordB2(e) >= bl[i][j].coordB &&
                                graph->getEdgeCoordE2(e) <= bl[i][j].coordE) {
                                return j;
                            }
                        }
                    }


                }
            }
            return -1;
        }
    }
}

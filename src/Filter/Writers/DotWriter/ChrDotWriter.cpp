#include <set>
#include <iomanip>
#include "ChrDotWriter.h"

namespace filter {
    namespace writers {
        void filter::writers::ChrDotWriter::writeVertexSet(std::vector<int> vert, std::string fileName) {
            TRACE("start write vert set");
            std::map<std::string, std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment> > > vertInChr;

            for (int v : vert) {
                std::vector<alig_info::InfoAboutContigsAlig::Alignment> aligs = aligInfo.getAlignment(v);
                for (auto alig : aligs) {
                    if ((alig.coordEnd - alig.coordBegin) * 10 > graph->getTargetLen(v)) {
                        vertInChr[alig.chrName].push_back(std::make_pair(v, alig));
                    }
                }
            }

            for (auto chrs : vertInChr) {
                if (chrs.first[chrs.first.size() - 1] != 'v') {
                    writeOneChr(chrs.first, chrs.second, chrs.first + "-rev", vertInChr[chrs.first + "-rev"], fileName);
                }
            }
            TRACE("finish write vert set");
        }

        bool cmp(std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > a,
                 std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > b) {
            return (a.second.coordBegin < b.second.coordBegin) ||
                    (a.second.coordBegin == b.second.coordBegin && a.second.chrName < b.second.chrName);
        }

        void ChrDotWriter::writeOneChr(const std::string chr,
                                       std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment >> verts,
                                       const std::string chr2,
                                       std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment >> verts2,
                                       std::string fileName) {

            TRACE("start write chr " << chr);

            fileName = fileName + "_" + chr;
            std::stable_sort(verts.begin(), verts.end(), cmp);
            std::stable_sort(verts2.rbegin(), verts2.rend(), cmp);

            assert(verts.size() == verts2.size());

            int pos = 0;
            int id = 0;
            while (pos < verts.size()) {
                std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > > chrv;
                std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > > chrv2;
                std::vector<int> toPrint;

                findVertsForCurFile(pos, verts, verts2, chrv, chrv2, toPrint);

                std::stringstream newFileName;
                newFileName << fileName << "_" << std::setfill('0') << std::setw(5) << id;

                writeOnePart(chr, chrv, chrv2, toPrint, newFileName.str());

                ++id;
            }

            TRACE("finish write chr " << chr);
        }

        void ChrDotWriter::writeOnePart(const std::string chrName,
                                        std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > > chrV,
                                        std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment > > chrV2,
                                        std::vector<int> allVert,
                                        std::string fileName) {
            TRACE("write vertex set");
            std::vector<bool> hasOtherEdge;
            std::vector<std::pair<int, int>> weightEdge;

            std::vector<int> vert;
            for (int v : allVert) {
                vert.push_back(v);
            }
            std::vector<int> coord;
            std::vector<int> coord2;
            std::vector<std::pair<int, int> >  vertE;
            std::vector<std::pair<int, int> > vertE2;

            for (auto v : chrV) {
                vert.push_back(v.first);
                coord.push_back(v.second.coordBegin / 1000);
                coord.push_back(std::max(v.second.coordEnd / 1000, v.second.coordBegin / 1000) );
                vertE.push_back(std::make_pair(std::max(v.second.coordEnd / 1000 , v.second.coordBegin / 1000), v.first));
            }

            for (auto v : chrV2) {
                vert.push_back(v.first);
                coord2.push_back(v.second.coordBegin / 1000);
                coord2.push_back(std::max(v.second.coordEnd / 1000, v.second.coordBegin / 1000) );
                vertE2.push_back(std::make_pair(std::max(v.second.coordEnd / 1000, v.second.coordBegin / 1000), v.first));
            }


            std::sort(coord.begin(), coord.end());
            coord.resize(std::unique(coord.begin(), coord.end()) - coord.begin());
            std::sort(vertE.begin(), vertE.end());

            std::sort(coord2.begin(), coord2.end());
            coord2.resize(std::unique(coord2.begin(), coord2.end()) - coord2.begin());
            std::reverse(coord2.begin(), coord2.end());
            std::sort(vertE2.rbegin(), vertE2.rend());

            hasOtherEdge.resize(vert.size());
            findVertWithOtherEdges(vert, hasOtherEdge, weightEdge);

            if (coord.size() == 0) return;

            std::ofstream out(fileName);

            out << "digraph {\n";
            out << "newrank=true;\n";
            out << "compound=true;\n";
            out << "outputorder=edgesfirst;\n";
            out << "ranksep=.75; size = \"7.5,7.5\";\n";

            writeOneStrand(chrV,
                           std::vector<bool>(hasOtherEdge.begin() + allVert.size(),
                                             hasOtherEdge.begin() + allVert.size() + chrV.size()), coord, vertE,
                           out);

            writeOneStrand(chrV2,
                           std::vector<bool>(hasOtherEdge.begin() + allVert.size() + chrV.size(),
                                             hasOtherEdge.end()), coord2, vertE2,
                           out);


            for (int i = 0; i < (int)std::min(coord.size(), coord2.size()); ++i) {
                out << "{ rank = same; " << coord[i] << "; " << coord2[i] << "; }\n";
            }

            for (int i = 0; i < (int)allVert.size(); ++i) {
                writeOneVertex(allVert[i], hasOtherEdge[i], out);
            }

            std::sort(weightEdge.rbegin(), weightEdge.rend());
            for (int i = 0; i < (int) weightEdge.size(); ++i) {
                writeOneEdge(weightEdge[i].second, out);
            }

            out << "labelloc=\"t\"" << "\n";
            out << "label=\"file name=" << fileName << "\";\n";

            out << "}\n";
            out.close();

            TRACE("finish write vert set");
        }

        std::ofstream & ChrDotWriter::writeOneStrand(std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment>> chrV,
                                                             std::vector<bool, std::allocator<bool>> hasOtherEdge, std::vector<int> coord,
                                                             std::vector<std::pair<int, int>> vertE, std::ofstream& out) {
            out << "subgraph {\n";

            out << "{\n";
            out << "node [shape=plaintext, fontsize=18];\n";
            out << coord[0];
            for (int i = 1; i < coord.size(); ++i) {
                out << " -> " << coord[i];
            }
            out << "\n";
            out << "}\n";

            out << "{\n";
            out << "node [shape=box];\n";

            for (int i = 0; i < (int) chrV.size(); ++i) {
                out << "subgraph cluster" << chrV[i].first << " {\n";
                out << "style=filled;\n";
                out << "color=\"#afdafc\";\n";

                this->writeOneVertex(chrV[i].first, hasOtherEdge[i], out);
                out << "e" << chrV[i].first << "[style=\"invis\", label=\"\", shape=point, fontsize=0, fixedsize=\"true\"];\n";
                out << "}\n";
            }

            int cur = 0, cur2 = 0;
            int ps = 0;
            while (cur < chrV.size() || cur2 < vertE.size()) {
                out << "{ rank = same; " << coord[ps] << "; ";
                while (cur < chrV.size() && chrV[cur].second.coordBegin/1000 == coord[ps]) {
                    out << "\"" << this->graph->getTargetName(chrV[cur].first) << "\"; ";
                    ++cur;
                }

                while (cur2 < vertE.size() && vertE[cur2].first == coord[ps]) {
                    out << "\"e" << vertE[cur2].second << "\"; ";
                    ++cur2;
                }

                ++ps;
                out << " }\n";
            }

            out << "}\n";



            out << "}\n";
            return out;
        }

        void ChrDotWriter::writeOneEdge(int e, std::ofstream &out) {
            TRACE("write one edge");

            int v = graph->getEdgeFrom(e);
            int u = graph->getEdgeTo(e);
            TRACE("write one edge e=" << e << " v=" << v << " u=" << u);

            out << "    \"" << graph->getTargetName(v) << "\" -> \"";
            out << graph->getTargetName(u) << "\" [ ";
            out << "color = \"" << graph->getLibColor(graph->getEdgeLib(e)) << "\", ";
            out << "penwidth = " << 1 + (int) log10(graph->getEdgeWeight(e)) << ", ";
            out << "label = " << "\"" << graph->getLibName(graph->getEdgeLib(e));
            out << "\n weight = " << (graph->getEdgeWeight(e));
            out << "\n id = " << e;
            out << "\n " << graph->getInfo(e) << "\" ];\n ";
            /*out << "ltail = cluster" << v << ", ";
            out << "lhead = cluster" << u << "];\n";*/
            //constraint=false
        }

        void ChrDotWriter::findVertsForCurFile(int &pos,
                                               std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment>> &verts,
                                               std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment>> &verts2,
                                               std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment>> &chrv,
                                               std::vector<std::pair<int, alig_info::InfoAboutContigsAlig::Alignment>> &chrv2,
                                               std::vector<int> &toPrint) {
            std::set<int> useVert;
            while (pos < verts.size() && useVert.size() < 5) {
                if (useVert.count(verts[pos].first) == 0) {
                    chrv.push_back(verts[pos]);
                    useVert.insert(verts[pos].first);
                    chrv2.push_back(verts2[pos]);
                    useVert.insert(verts2[pos].first);
                }
                ++pos;
            }

            for (auto v : chrv) {
                std::vector<int> nw = searcher.findVertInLocalArea(v.first, 1);
                if (nw.size() > 5) continue;
                for (int u : nw) {
                    if (useVert.count(u) == 0) {
                        toPrint.push_back(u);
                    }
                }
            }


            for (auto v : chrv2) {
                std::vector<int> nw = searcher.findVertInLocalArea(v.first, 1);
                if (nw.size() > 5) continue;
                for (int u : nw) {
                    if (useVert.count(u) == 0) {
                        toPrint.push_back(u);
                    }
                }
            }
        }
    }
}

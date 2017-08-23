#include <set>
#include <iomanip>
#include "ChrDotWriter.h"

namespace filter {
    namespace writers {
        void filter::writers::ChrDotWriter::writeVertexSet(std::vector<int> vert, std::string fileName) {
            TRACE("start write vert set");
            std::map<std::string, std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment> > > vertInChr;

            for (int v : vert) {
                std::vector<statistics::InfoAboutContigsAlig::Alignment> aligs = aligInfo.getAlignment(v);
                for (auto alig : aligs) {
                    if ((alig.coordEnd - alig.coordBegin) * 10 > graph->getTargetLen(v)) {
                        vertInChr[alig.chrName].push_back(std::make_pair(v, alig));
                    }
                }
            }

            for (auto chrs : vertInChr) {
                writeOneChr(chrs.first, chrs.second, fileName);
            }
            TRACE("finish write vert set");
        }

        bool cmp(std::pair<int, statistics::InfoAboutContigsAlig::Alignment > a,
                 std::pair<int, statistics::InfoAboutContigsAlig::Alignment > b) {
            return a.second.coordBegin < b.second.coordBegin;
        }

        void ChrDotWriter::writeOneChr(const std::string chr,
                                       std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment >> verts,
                                       std::string fileName) {
            TRACE("start write chr " << chr);
            fileName = fileName + "_" + chr;
            int step = 20;
            int cur = 0;

            std::sort(verts.begin(), verts.end(), cmp);

            int pos = 0;
            int id = 0;
            while (pos < verts.size()) {
                std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment > > chrv;
                std::set<int> useVert;

                while (pos < verts.size() && useVert.size() < 20) {
                    if (useVert.count(verts[pos].first) == 0) {
                        chrv.push_back(verts[pos]);
                        useVert.insert(verts[pos].first);
                    }
                    ++pos;
                }

                std::vector<int> toPrint;
                for (auto v : chrv) {
                    std::vector<int> nw = searcher.findVertInLocalArea(v.first, 1);
                    if (nw.size() > 20) continue;
                    for (int u : nw) {
                        if (useVert.count(u) == 0) {
                            toPrint.push_back(u);
                        }
                    }
                }

                std::stringstream newFileName;
                newFileName << fileName << "_" << std::setfill('0') << std::setw(5) << id;

                writeOnePart(chr, chrv, toPrint, newFileName.str());

                ++id;
            }

            TRACE("finish write chr " << chr);
        }

        void ChrDotWriter::writeOnePart(const std::string chrName, std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment > > chrV, std::vector<int> allVert,
                                        std::string fileName) {
            TRACE("write vertex set");
            std::vector<bool> hasOtherEdge;
            std::vector<std::pair<int, int>> weightEdge;

            std::vector<int> vert;
            for (int v : allVert) {
                vert.push_back(v);
            }
            std::vector<int> coord;

            for (auto v : chrV) {
                vert.push_back(v.first);
                if (coord.size() == 0 || coord[coord.size() - 1] != v.second.coordBegin/1000) {
                    coord.push_back(v.second.coordBegin / 1000);
                }
            }

            hasOtherEdge.resize(vert.size());
            findVertWithOtherEdges(vert, hasOtherEdge, weightEdge);

            if (coord.size() == 0) return;

            std::ofstream out(fileName);

            out << "digraph {\n";
            out << "ranksep=.75; size = \"7.5,7.5\";\n";
            out << "{\n";
            out << "node [shape=plaintext, fontsize=16];\n";
            out << coord[0];
            for (int i = 1; i < coord.size(); ++i) {
                out << " -> " << coord[i];
            }
            out << "\n";
            out << "}\n";

            out << "{\n";
            out << "node [shape=box];\n";

            for (int i = 0; i < (int) chrV.size(); ++i) {
                writeOneVertex(chrV[i].first, hasOtherEdge[allVert.size() + i], out);
            }

            int cur = 0;
            int ps = 0;
            while (cur < chrV.size()) {
                out << "{ rank = same; " << coord[ps] << "; ";
                TRACE(ps << " " << chrV[cur].second.coordBegin/1000);
                while (cur < chrV.size() && chrV[cur].second.coordBegin/1000 == coord[ps]) {
                    out << "\"" << graph->getTargetName(chrV[cur].first) << "\"; ";
                    ++cur;
                }

                ++ps;
                out << " }\n";
            }

            out << "}\n";

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
    }
}

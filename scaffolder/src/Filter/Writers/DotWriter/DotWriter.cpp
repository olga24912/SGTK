#include <iostream>
#include "DotWriter.h"
#include <cmath>

namespace filter {
    namespace writers {
        void DotWriter::writeVertexSet(std::vector<int> vert, std::string fileName) {
            DEBUG("writeVertexSet");
            std::vector<std::vector<int> > res = graphSplitter.split(filter, vert);
            for (int i = 0; i < (int) res.size(); ++i) {
                std::stringstream ss;
                ss << fileName << i;
                std::string name = std::string(ss.str());
                if (validator->isGoodVertexSet(res[i], filter)) {
                    TRACE("isOk vertex set" << i)
                    writeOneVertexSet(res[i], name);
                }
            }
        }

        void DotWriter::writeOneVertex(int v, bool isColored, std::ofstream &out) {
            TRACE("write one vertex v=" << v << " isColored=" << isColored);

            out << "    \"" << filter->getTargetName(v) << "\"[label=\" " << filter->getTargetName(v) <<
                " id = " << v
                << "\nlen = " << filter->getTargetLen(v) << "\"";
            if (isColored) {
                out << " , style = \"filled\", color = \"#F0E68C\"";
            }
            out << "];\n";
        }

        void DotWriter::writeOneEdge(int e, std::ofstream &out) {
            TRACE("write one edge e=" << e);

            int v = filter->getEdgeFrom(e);
            int u = filter->getEdgeTo(e);
            out << "    \"" << filter->getTargetName(v) << "\" -> \"";
            out << filter->getTargetName(u) << "\" [ ";
            out << "color = \"" << filter->getLibColor(filter->getEdgeLib(e)) << "\", ";
            out << "penwidth = " << 1 + (int) log10(filter->getEdgeWeight(e)) << ", ";
            out << "label = " << "\"" << filter->getLibName(filter->getEdgeLib(e));
            out << "\n weight = " << (filter->getEdgeWeight(e));
            out << "\n id = " << e;
            out << "\n " << filter->getInfo(e) << "\" ]\n";
        }

        void DotWriter::writeOneVertexSet(std::vector<int> vert, std::string fileName) {
            TRACE("write vertex set");
            std::vector<bool> hasOtherEdge((unsigned) filter->getVertexCount(), 0);
            std::vector<std::pair<int, int> > weightEdge;
            for (int i = 0; i < (int) vert.size(); ++i) {
                int v = vert[i];
                for (int e : filter->getEdges(v)) {
                    int u = filter->getEdgeTo(e);
                    int was = 0;
                    for (int h = 0; h < (int) vert.size(); ++h) {
                        if (vert[h] == u) was = 1;
                    }
                    if (was) {
                        weightEdge.push_back(std::make_pair(filter->getEdgeWeight(e), e));
                    } else {
                        hasOtherEdge[i] = 1;
                    }
                }
                for (int e : filter->getEdgesR(v)) {
                    int u = filter->getEdgeFrom(e);
                    int was = 0;
                    for (int h = 0; h < (int) vert.size(); ++h) {
                        if (vert[h] == u) was = 1;
                    }
                    if (!was) {
                        hasOtherEdge[i] = 1;
                    }
                }
            }
            if (vert.size() == 1) return;

            std::ofstream out(fileName);
            out << "digraph {\n";
            for (int i = 0; i < (int) vert.size(); ++i) {
                writeOneVertex(vert[i], hasOtherEdge[i], out);
            }
            std::sort(weightEdge.rbegin(), weightEdge.rend());
            for (int i = 0; i < (int) weightEdge.size(); ++i) {
                writeOneEdge(weightEdge[i].second, out);
            }
            out << "}\n";
            out.close();
        }
    }
}
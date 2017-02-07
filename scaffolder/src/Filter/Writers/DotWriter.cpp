#include <iostream>
#include "DotWriter.h"

void DotWriter::writeVertexSet(std::vector<int> vert, std::string fileName) {
    std::cerr << "write set: "  << vert.size() << std::endl;
    std::vector<std::vector<int> > res = graphSplitter.split(filter, vert);
    for (int i = 0; i < (int)res.size(); ++i) {
        std::stringstream ss;
        ss << fileName << i;
        std::string name = std::string(ss.str());
        writeOneVertexSet(res[i], name);
    }
}

void DotWriter::writeOneVertex(int v, bool isColored, std::ofstream &out) {
    out << "    \"" << filter->getTargetName(v) << "\"[label=\" " << filter->getTargetName(v) <<
        " id = " << v
        << "\nlen = " << filter->getTargetLen(v) << "\"";
    if (isColored) {
        out << " , style = \"filled\", color = \"#F0E68C\"";
    }
    out << "];\n";
}

void DotWriter::writeOneEdge(int e, std::ofstream &out) {
    int v = filter->getEdgeFrom(e);
    int u = filter->getEdgeTo(e);
    out << "    \"" << filter->getTargetName(v) << "\" -> \"";
    out << filter->getTargetName(u) << "\" [ ";
    out << "color = \"" << filter->getLibColor(filter->getEdgeLib(e)) << "\", ";
    out << "penwidth = "<< 1 + (int)log10(filter->getEdgeWieght(e)) << ", ";
    out << "label = " << "\"" << filter->getLibName(filter->getEdgeLib(e));
    out << "\n weight = " << (filter->getEdgeWieght(e));
    out << "\n id = "<< e << "\" ]\n";
}

void DotWriter::writeOneVertexSet(std::vector<int> vert, std::string fileName) {
    std::vector<bool> hasOtherEdge((unsigned)filter->getVertexCount(), 0);
    std::vector<std::pair<int, int> > weightEdge;
    for (int i = 0; i < (int)vert.size(); ++i) {
        int v = vert[i];
        for (int e : filter->getEdges(v)) {
            int u = filter->getEdgeTo(e);
            int was = 0;
            for (int h = 0; h < (int)vert.size(); ++h) {
                if (vert[h] == u) was = 1;
            }
            if (was) {
                weightEdge.push_back(std::make_pair(filter->getEdgeWieght(e), e));
            } else {
                hasOtherEdge[i] = 1;
            }
        }
        for (int e : filter->getEdgesR(v)) {
            int u = filter->getEdgeFrom(e);
            int was = 0;
            for (int h = 0; h < (int)vert.size(); ++h) {
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
    for (int i = 0; i < (int)vert.size(); ++i) {
        writeOneVertex(vert[i], hasOtherEdge[i], out);
    }
    std::sort(weightEdge.rbegin(), weightEdge.rend());
    for(int i = 0; i < (int)weightEdge.size(); ++i) {
        writeOneEdge(weightEdge[i].second, out);
    }
    out << "}\n";
    out.close();
}

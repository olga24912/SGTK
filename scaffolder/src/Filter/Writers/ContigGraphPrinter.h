//
// Created by olga on 23.10.16.
//

/*#ifndef SCAFFOLDER_CONTIGGRAPHPRINTER_H
#define SCAFFOLDER_CONTIGGRAPHPRINTER_H

#include "ContigGraph/ContigGraph.h"

using namespace std;

class ContigGraphPrinter {
private:
    struct IsGoodEdge {
    public:
        IsGoodEdge(){
            libNum = -1;
        }
        IsGoodEdge(int libNum) {
            this->libNum = libNum;
        }
        bool operator()(ContigGraph *g, int e) const {
            if (g->edgeWeight[e] < g->libMinEdgeWight[g->edgeLib[e]]) {
                return false;
            }
            if (libNum != -1) {
                if (g->edgeLib[e] != libNum) return false;
            }
            return true;
        }

    private:
        int libNum;
    };

    static vector<int> vertexInBigComponents(ContigGraph *g, int size);

public:
    static void writeFullGraphDotFormat(ContigGraph * graph, string fileName);
    static void writeAllLocalGraphDotFormat(ContigGraph *graph, int dist);
    static void writeLocalGraph(ContigGraph *g, int dist, int v, string fileName);
    static void writeLocalSegGraph(ContigGraph *g, int dist, int vb, int ve, string fileName);
    static void writeBigComponent(ContigGraph *g, int minSize, string fileName);
    static void writeSplitBigComponent(ContigGraph *g, int minSize, string fileName);
    static void writeAlongPath(ContigGraph *g, int libId, int dist, int minSize, string fileName);
};


#endif //SCAFFOLDER_CONTIGGRAPHPRINTER_H
*/
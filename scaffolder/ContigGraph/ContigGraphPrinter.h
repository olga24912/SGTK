//
// Created by olga on 23.10.16.
//

#ifndef SCAFFOLDER_CONTIGGRAPHPRINTER_H
#define SCAFFOLDER_CONTIGGRAPHPRINTER_H

#include "ContigGraph.h"

using namespace std;

class ContigGraphPrinter {
private:
    static void writeOneEdge(ContigGraph * graph, ofstream& out, int v, int e);
    static void writeOneVertex(ContigGraph * g, ofstream& out, int v);
    static vector<int> findAllVert(ContigGraph *g, int dist, int v);
    static void dfsFindComponent(ContigGraph *g, int* color, int currentCol, int v);
    static vector<int> vertexInBigComponents(ContigGraph *g, int size);
    static void writeThisVertex(ContigGraph *g, vector<int> &drawV, string fileName);

public:
    static void writeFullGraphDotFormat(ContigGraph * graph, string fileName);
    static void writeAllLocalGraphDotFormat(ContigGraph *graph, int dist);
    static void writeLocalGraph(ContigGraph *g, int dist, int v, string fileName);
    static void writeBigComponent(ContigGraph *g, int minSize, string fileName);
};


#endif //SCAFFOLDER_CONTIGGRAPHPRINTER_H

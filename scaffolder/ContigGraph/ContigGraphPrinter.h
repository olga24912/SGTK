//
// Created by olga on 23.10.16.
//

#ifndef SCAFFOLDER_CONTIGGRAPHPRINTER_H
#define SCAFFOLDER_CONTIGGRAPHPRINTER_H

#include "ConigGraph.h"

using namespace std;

class ContigGraphPrinter {
    static void writeOneEdge(ConigGraph * graph, ofstream& out, int v, int e);
    static void writeOneVertex(ConigGraph * g, ofstream& out, int v);
    static vector<int> findAllVert(ConigGraph *g, int dist, int v);

public:
    static void writeFullGraphDotFormat(ConigGraph * graph, string fileName);
    static void writeAllLocalGraphDotFormat(ConigGraph *graph, int dist);

    static void writeLocalGraph(ConigGraph *g, int dist, int v, string fileName);
};


#endif //SCAFFOLDER_CONTIGGRAPHPRINTER_H

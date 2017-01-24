//
// Created by olga on 24.01.17.
//

#ifndef SCAFFOLDER_GRAPHWITHFILTER_H
#define SCAFFOLDER_GRAPHWITHFILTER_H

#include "ContigGraph/ContigGraph.h"

class GraphWithFilter {
    friend class ContigGraphPrinter;
private:
    ContigGraph *graph;
    int minContigLen;

    vector<bool> ignore;
    vector<int> libMinEdgeWight;
public:
    void filterByEdgeWeight(int minEdgeWeight); //delete edge with small wight

    void filterByContigLen(int minContigLen);//fileter vertex with small len

    void setIgnore(int v) {
        ignore[v] = ignore[v] ^ 1;
    }

    bool isIgnore(int v) {
        return ignore[v];
    }
};


#endif //SCAFFOLDER_GRAPHWITHFILTER_H

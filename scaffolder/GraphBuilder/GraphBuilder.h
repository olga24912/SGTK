//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H


#include "../ContigGraph/ContigGraph.h"

/*
 * main class for generate conection between contigs.
 */
class GraphBuilder {
protected:
    ContigGraph* graph;

    string libName;
    int minContigLen = 0;
    int minEdgeWight = 0;
public:
    /*
     * fun that need to call for add conection
     * between contigs;
     */
    virtual void evaluate() = 0;

    /*
     * set barrier contig len,
     * Contigs with smaller len will be ignore.
     */
    void setMinContigLen(int minContigLen);

    /*
     * set barrier edge wight.
     * Edge with smaller wieght will be remove.
     */
    void setMinEdgeWight(int minEdgeWight);
    void setGraph(ContigGraph* graph);

    void setLibName(string libName);
};


#endif //SCAFFOLDER_GRAPHBUILDER_H

//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H


#include "../ConigGraph.h"


class GraphBuilder {
protected:
    ConigGraph* graph;

    int minContigLen = 0;
    int minEdgeWight = 0;
public:
    virtual void evaluate() = 0;

    void setMinContigLen(int minContigLen);
    void setMinEdgeWight(int minEdgeWight);
    void setGraph(ConigGraph* graph);
};


#endif //SCAFFOLDER_GRAPHBUILDER_H

//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H


#include "../Graph.h"

class GraphBuilder {
protected:
    Graph graph;

    int minContigLen = 0;
    int minEdgeWight = 0;
public:
    virtual void evaluate() = 0;

    void setMinContigLen(int minContigLen);
    void setMinEdgeWight(int minEdgeWight);
    void setGraph(Graph &graph);
};


#endif //SCAFFOLDER_GRAPHBUILDER_H

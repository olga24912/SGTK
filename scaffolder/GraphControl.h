//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHCONTROL_H
#define SCAFFOLDER_GRAPHCONTROL_H

#include "Graph.h"
#include "GraphBuilder/DNAPairReadGraphBuilder.h"

class GraphControl {
private:
    Graph graph;
public:
    void evaluate(int argc, char **argv);
};


#endif //SCAFFOLDER_GRAPHCONTROL_H

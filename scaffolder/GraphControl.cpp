//
// Created by olga on 08.10.16.
//

#include "GraphControl.h"

void GraphControl::evaluate(int argc, char **argv) {
    DNAPairReadGraphBuilder gb;
    int pos = 1;
    while (pos < argc) {
        cerr << pos << " " << argc << endl;

        gb.setMinContigLen(atoi(argv[pos + 4]));
        gb.setFileName1(argv[pos]);
        gb.setFileName2(argv[pos + 1]);
        gb.setDistBetweenPairReads(atoi(argv[pos + 3]));
        gb.setMinEdgeWight(atoi(argv[pos + 2]));


        graph.newLib();
        gb.setGraph(&graph);
        gb.evaluate();

        pos += 5;
    }

    graph.writeGraphDotFormat("graph.dot");
}

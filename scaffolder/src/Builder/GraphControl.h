#ifndef SCAFFOLDER_GRAPHCONTROL_H
#define SCAFFOLDER_GRAPHCONTROL_H

#include "ContigGraph/ContigGraph.h"
#include "Builder/GraphBuilder/DNAPairReadGraphBuilder.h"
#include "Builder/GraphBuilder/PairReadGraphBuilder.h"

namespace builder {
// that class parse args and call other class for
// build conection and keep info about conection.
    class GraphControl {
    private:
        contig_graph::ContigGraph graph;
        std::string path = ".";
    public:
        void evaluate(int argc, char **argv); //build conection with this args

    private:
        DECL_LOGGER("GraphControl");
    };
}

#endif //SCAFFOLDER_GRAPHCONTROL_H

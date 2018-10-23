#ifndef SCAFFOLDER_GRAPHCONTROL_H
#define SCAFFOLDER_GRAPHCONTROL_H

#include "ContigGraph/ContigGraph.h"
#include "Builder/GraphBuilder/PairReadGraphBuilder.h"

namespace builder {

/* Graph Control parse args and call other class for
* build connections and keep info about connections.
*/
    class GraphControl {
    private:
        contig_graph::ContigGraph graph;
        std::string path = ".";
    public:
        //build connection with this args
        void evaluate(int argc, char **argv);

    private:
        DECL_LOGGER("GraphControl");
    };
}

#endif //SCAFFOLDER_GRAPHCONTROL_H

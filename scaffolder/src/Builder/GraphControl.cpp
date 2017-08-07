#include "GraphControl.h"
#include "Builder/GraphBuilder/ReferenceGraphBuilder.h"

namespace builder {
    void GraphControl::evaluate(int argc, char **argv) {
        using namespace graph_builder;
        using namespace contig_graph;
        INFO("start build graph");

        if (argv[1] == "RNA_PAIR") {
            RNAPairReadGraphBuilder gb = RNAPairReadGraphBuilder();

            gb.setFileName1(argv[2]);
            gb.setFileName2(argv[3]);

            gb.setGraph(&graph);
            gb.evaluate();
        } else if (argv[1] == "DNA_PAIR") {
            DNAPairReadGraphBuilder gb = DNAPairReadGraphBuilder();

            gb.setFileName1(argv[2]);
            gb.setFileName2(argv[3]);

            gb.setDistBetweenPairReads(atoi(argv[4]));

            gb.setGraph(&graph);
            gb.evaluate();
        } else if (argv[1] == "RNA_SPLIT") {
            RNAPairReadGraphBuilder gb = RNAPairReadGraphBuilder();

            gb.setFileName1(argv[2]);
            gb.setFileName2(argv[3]);

            gb.setOneSideReadFlag(true);

            gb.setGraph(&graph);
            gb.evaluate();
        } else if (argv[1] == "REF") {
            ReferenceGraphBuilder gb = ReferenceGraphBuilder();

            gb.setRefFileName(argv[2]);
            gb.setQueryFileName(argv[3]);
            gb.setMinContigLen(atoi(argv[4]));

            gb.setGraph(&graph);
            gb.evaluate();
        } else {
            ERROR("unkonwn connection type: " << argv[1]);
        }

        graph.write(path + "/graph.gr");
        INFO("end build graph");
    }
}
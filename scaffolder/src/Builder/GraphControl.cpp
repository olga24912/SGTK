#include <Builder/GraphBuilder/RNAPairReadGraphBuilder.h>
#include <Builder/GraphBuilder/RNASplit50GraphBuilder.h>
#include <Builder/GraphBuilder/RNASplit30GraphBuilder.h>
#include "GraphControl.h"
#include "Builder/GraphBuilder/ReferenceGraphBuilder.h"

namespace builder {
    void GraphControl::evaluate(int argc, char **argv) {
        using namespace graph_builder;
        using namespace contig_graph;
        INFO("start build graph with type: " << argv[1]);

        if (std::string(argv[1]) == "RNA_PAIR") {
            RNAPairReadGraphBuilder* gb = new RNAPairReadGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);
            gb->setLibName(argv[4]);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "DNA_PAIR") {
            DNAPairReadGraphBuilder* gb = new DNAPairReadGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setDistBetweenPairReads(atoi(argv[4]));
            gb->setLibName(argv[5]);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "RNA_SPLIT_50") {
            RNASplit50GraphBuilder* gb = new RNASplit50GraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setDirect(atoi(argv[5]));

            gb->setGraph(&graph);
            gb->evaluate();
        } else if ( std::string(argv[1]) == "RNA_SPLIT_30") {
            RNASplit30GraphBuilder* gb = new RNASplit30GraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setDirect(atoi(argv[5]));

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "REF") {
            ReferenceGraphBuilder* gb = new ReferenceGraphBuilder;

            gb->setRefFileName(argv[2]);
            gb->setQueryFileName(argv[3]);
            gb->setMinContigLen(atoi(argv[4]));
            gb->setLibName(argv[5]);

            gb->setGraph(&graph);
            gb->evaluate();
        } else {
            ERROR("unkonwn connection type: " << argv[1]);
        }

        graph.write(path + "/graph.gr");
        INFO("end build graph");
    }
}
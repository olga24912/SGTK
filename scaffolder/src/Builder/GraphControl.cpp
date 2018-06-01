#include <Builder/GraphBuilder/RNASplitGraphBuilder.h>
#include <Builder/GraphBuilder/ConnectionGraphBuilder.h>
#include <Builder/GraphBuilder/DNAMatePairReadsGraphBuilder.h>
#include <Builder/GraphBuilder/PacbioGraphBuilder.h>
#include <Builder/GraphBuilder/FASTGGraphBuilder.h>
#include "GraphControl.h"
#include "Builder/GraphBuilder/ReferenceGraphBuilder.h"

namespace builder {
    void GraphControl::evaluate(int argc, char **argv) {
        using namespace graph_builder;
        using namespace contig_graph;

        if (std::string(argv[1]) == "RNA_PAIR") {
            PairReadGraphBuilder* gb = new PairReadGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);
            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::RNA_PAIR);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "DNA_PAIR_FR") {
            DNAMatePairReadsGraphBuilder* gb = new DNAMatePairReadsGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::DNA_PAIR);

            gb->setRevFirstFlag(0);
            gb->setRevSecondFlag(1);
            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "DNA_PAIR_RF") {
            DNAMatePairReadsGraphBuilder* gb = new DNAMatePairReadsGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::DNA_PAIR);

            gb->setRevFirstFlag(1);
            gb->setRevSecondFlag(0);
            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "DNA_PAIR_FF") {
            DNAMatePairReadsGraphBuilder* gb = new DNAMatePairReadsGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::DNA_PAIR);

            gb->setRevFirstFlag(0);
            gb->setRevSecondFlag(0);
            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "PACBIO") {
            PacbioGraphBuilder* gb = new PacbioGraphBuilder;

            gb->setFileName(argv[2]);
            gb->setContigFileName(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::PACBIO);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "RNA_SPLIT_50") {
            RNASplitGraphBuilder* gb = new RNASplitGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::RNA_SPLIT_50);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if ( std::string(argv[1]) == "RNA_SPLIT_30") {
            RNASplitGraphBuilder* gb = new RNASplitGraphBuilder;

            gb->setFileName1(argv[2]);
            gb->setFileName2(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::RNA_SPLIT_30);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "REF") {
            ReferenceGraphBuilder* gb = new ReferenceGraphBuilder;

            gb->setRefFileName(argv[2]);
            gb->setQueryFileName(argv[3]);
            gb->setMinContigLen(atoi(argv[4]));
            gb->setLibName(argv[5]);
            gb->setLibType(ContigGraph::Lib::Type::REF);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "CONNECTION") {
            ConnectionGraphBuilder* gb = new ConnectionGraphBuilder;

            gb->setConnectionFile(argv[2]);
            gb->setContigFile(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::CONNECTION);

            gb->setGraph(&graph);
            gb->evaluate();
        } else if (std::string(argv[1]) == "FASTG") {
            FASTGGraphBuilder* gb = new FASTGGraphBuilder;

            gb->setFASTGFile(argv[2]);
            gb->setContigFile(argv[3]);

            gb->setLibName(argv[4]);
            gb->setLibType(ContigGraph::Lib::Type::FASTG);

            gb->setGraph(&graph);
            gb->evaluate();
        } else {
            ERROR("unkonwn connection type: " << argv[1]);
        }

        graph.write(path + "/graph.gr");
    }
}
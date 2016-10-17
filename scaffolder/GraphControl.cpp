//
// Created by olga on 08.10.16.
//

#include "GraphControl.h"
#include "GraphBuilder/RNAPairReadGraphBuilder.h"
#include "GraphBuilder/RNASplitReadGraphBuilder.h"

void GraphControl::evaluate(int argc, char **argv) {
    argc -= (argc > 0);
    argv += (argc > 0);

    option::Stats stats(usage, argc, argv);
    option::Option options[stats.options_max], buffer[stats.buffer_max];

    option::Parser parse(usage, argc, argv, options, buffer);

    if (parse.error() || options[HELP] || argc == 0) {
        option::printUsage(std::cout, usage);
        return;
    }

    GraphBuilder* gb = nullptr;
    bool needEval = false;

    for (int i = 0; i < parse.optionsCount(); ++i) {
        option::Option& opt = buffer[i];
        string arg;
        switch (opt.index()) {
            case NEW:
                if (needEval) {
                    graph.newLib();
                    gb -> setGraph(&graph);
                    gb -> evaluate();
                }

                delete gb;
                needEval = true;

                arg = opt.arg;
                if (arg == "RNA_PAIR") {
                    gb = new RNAPairReadGraphBuilder();
                } else if (arg == "DNA_PAIR") {
                    gb = new DNAPairReadGraphBuilder();
                } else if (arg == "RNA_SPLIT") {
                    gb = new RNASplitReadGraphBuilder();
                }

                break;
            case SAM1:
                if ( dynamic_cast <PairReadGraphBuilder *> ( gb )) {
                    (dynamic_cast<PairReadGraphBuilder*> (gb)) -> setFileName1(opt.arg);
                } else {
                    printf("Put the pair read sam file for not pair read builder.");
                    return;
                }
                break;
            case SAM2:
                if ( dynamic_cast <PairReadGraphBuilder *> ( gb )) {
                    (dynamic_cast<PairReadGraphBuilder*> (gb)) -> setFileName2(opt.arg);
                } else {
                    printf("Put the pair read sam file for not pair read builder.");
                    return;
                }
                break;
            case MINCONTIGLEN:
                gb->setMinContigLen(atoi(opt.arg));
                break;
            case MINEDGELEN:
                gb->setMinEdgeWight(atoi(opt.arg));
                break;
            default:break;
        }
    }

    if (needEval) {
        graph.newLib();
        gb -> setGraph(&graph);
        cerr << "ok";
        gb -> evaluate();
    }

    delete gb;

/*
    int pos = 1;
    while (pos < argc) {
        //cerr << pos << " " << argc << endl;

        //gb.setMinContigLen(atoi(argv[pos + 3]));
        //gb.setFileName1(argv[pos]);
        //gb.setFileName2(argv[pos + 1]);
        //gb.setDistBetweenPairReads(atoi(argv[pos + 3]));
        //gb.setMinEdgeWight(atoi(argv[pos + 2]));

        //gb.setRefFileName(argv[pos]);
        //gb.setRnaReadFileName(argv[pos + 1]);

        //*graph.newLib();
        //gb.setGraph(&graph);
        //gb.evaluate();

        pos += 2;
    }*/

    graph.writeGraphDotFormat("graph.dot");
}

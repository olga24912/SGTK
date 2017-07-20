#include <algorithm>
#include <iostream>
#include <set>
#include "ScaffolderPipeline.h"
#include "ScaffoldStrategy.h"
#include "ScaffoldStrategyUniqueConnection.h"
#include "RuleBigDifInWeight.h"
#include "RuleDelCycle.h"
#include "RuleInOneLine.h"
#include "RuleBigDeg.h"
#include "RuleDelSmallCycle.h"

namespace filter {
    namespace scaffolder {
        void ScaffolderPipeline::evaluate(ContigGraph *graph, std::string contigFile, std::string out) {
            INFO("start build scaffolds");
            std::vector<int> libs = graph->getLibList();
            for (int i = (int)libs.size() - 1; i > 0; --i) {
                double w1 = 1, w2 = 1;
                if (i == libs.size() - 1 && graph->getLibType(libs[i]) == contig_graph::ContigGraph::Lib::RNA_SPLIT_50) {
                    w1 = 1.75;
                }
                if (graph->getLibType(libs[i - 1]) == contig_graph::ContigGraph::Lib::RNA_SPLIT_50) {
                    w2 = 1.75;
                }
                graph->mergeLib(libs[i], libs[i - 1], "lib", w1, w2);
            }

            Scaffolds scaffolds(contigFile);

            RuleDelSmallCycle rdsc;
            rdsc.simplifyGraph(graph);
            RuleBigDifInWeight rbd;
            rbd.simplifyGraph(graph);
            RuleBigDeg rbdeg;
            rbdeg.simplifyGraph(graph);
            RuleInOneLine riol;
            riol.simplifyGraph(graph);
            RuleDelCycle rdc;
            rdc.simplifyGraph(graph);

            graph->write("smp.gr");

            ScaffoldStrategyUniqueConnection ssuc;
            ssuc.addConnection(&scaffolds, graph);
            scaffolds.print(out);
        }
    }
}

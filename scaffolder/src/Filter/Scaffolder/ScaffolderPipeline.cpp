#include <algorithm>
#include <iostream>
#include <set>
#include <Filter/CommandParsers/Commands/Command.h>
#include "ScaffolderPipeline.h"
#include "ScaffoldStrategy.h"
#include "ScaffoldStrategyUniqueConnection.h"
#include "RuleBigDifInWeight.h"
#include "RuleDelCycle.h"
#include "RuleInOneLine.h"
#include "RuleBigDeg.h"
#include "RuleDelSmallCycle.h"
#include "RuleDelSmallEdges.h"
#include "RuleCoord.h"
#include "RuleDel30.h"
#include "RuleValidateCoord.h"
#include "RuleCovering.h"
#include "RuleExonBlocks.h"

namespace filter {
    namespace scaffolder {
        using namespace commands;
        void ScaffolderPipeline::evaluate(ContigGraph *graph, std::string contigFile,
                                          std::string out, std::vector<State::BamFiles> bamFiles) {
            INFO("start build scaffolds");

            Scaffolds scaffolds(contigFile);

            RuleDelSmallEdges rdse;
            rdse.simplifyGraph(graph);
            graph->write("smp.gr");
            RuleDelSmallCycle rdsc;
            rdsc.simplifyGraph(graph);
            RuleInOneLine riol;
            riol.simplifyGraph(graph);
            RuleBigDifInWeight rbddiw;
            rbddiw.simplifyGraph(graph);
            riol.simplifyGraph(graph);
            RuleExonBlocks reb;
            reb.simplifyGraph(graph);
            RuleDel30 rd30;
            rd30.simplifyGraph(graph);

            ScaffoldStrategyUniqueConnection ssuc;
            ssuc.addConnection(&scaffolds, graph);
            scaffolds.print(out);
        }
    }
}

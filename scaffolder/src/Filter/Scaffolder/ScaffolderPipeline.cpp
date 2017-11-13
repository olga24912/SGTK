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
#include "ScaffoldStrategyOneLine.h"

namespace filter {
    namespace scaffolder {
        using namespace commands;
        void ScaffolderPipeline::evaluate(ContigGraph *graph, std::string contigFile,
                                          std::string out, std::vector<State::BamFiles> bamFiles) {
            Scaffolds scaffolds(contigFile);


            INFO("Simplification step 1: delete small edges");
            RuleDelSmallEdges rdse;
            rdse.simplifyGraph(graph);

            INFO("Simplification step 2: delete cycles");
            RuleDelSmallCycle rdsc;
            rdsc.simplifyGraph(graph);

            INFO("Simplification step 3: fork with big diff in wieght");
            RuleBigDifInWeight rbddiw;
            rbddiw.simplifyGraph(graph);

            INFO("Simplification step 4: validation by gene annotation");
            RuleExonBlocks reb;
            reb.simplifyGraph(graph);

            INFO("Simplification step 5: delete split-30 edges");
            RuleDel30 rd30;
            rd30.simplifyGraph(graph);

            INFO("Build scaffolds");
            ScaffoldStrategyOneLine ssuc;
            ssuc.addConnection(&scaffolds, graph);

            INFO("Print scaffold to " + out);
            scaffolds.print(out);
        }
    }
}

#ifndef SCAFFOLDER_RULEDELSMALLCYCLE_H
#define SCAFFOLDER_RULEDELSMALLCYCLE_H

#include <ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleDelSmallCycle : public Rule {
        private:
            const int MAX_DIF = 5;
        public:
            void simplifyGraph(ContigGraph *graph) override;

        };
    }
}

#endif //SCAFFOLDER_RULEDELSMALLCYCLE_H

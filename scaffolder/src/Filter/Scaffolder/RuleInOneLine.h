#ifndef SCAFFOLDER_RULEINONELINE_H
#define SCAFFOLDER_RULEINONELINE_H

#include <ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleInOneLine : Rule {
        public:
            void simplifyGraph(ContigGraph *graph) override;

        private:
            bool havePath(ContigGraph *filter, int u, int w);
        };
    }
}


#endif //SCAFFOLDER_RULEINONELINE_H

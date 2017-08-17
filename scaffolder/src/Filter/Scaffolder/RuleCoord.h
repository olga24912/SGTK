#ifndef SCAFFOLDER_RULECOORD_H
#define SCAFFOLDER_RULECOORD_H

#include <Filter/ContigGraph/ContigGraph.h>
#include <set>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;

        class RuleCoord : public Rule {
        private:
            std::set<int> edgeForDel;
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void directSimpl(ContigGraph *graph, int v);

            void revSimpl(ContigGraph *graph, int v);

            int cntEdge(ContigGraph *graph, std::vector<int> edges);
        };
    }
}



#endif //SCAFFOLDER_RULECOORD_H

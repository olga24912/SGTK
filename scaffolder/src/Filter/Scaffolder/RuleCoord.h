#ifndef SCAFFOLDER_RULECOORD_H
#define SCAFFOLDER_RULECOORD_H

#include <Filter/ContigGraph/ContigGraph.h>
#include <set>
#include <Filter/Statistics/InfoAboutContigsAlig.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;

        class RuleCoord : public Rule {
        private:
            std::set<int> edgeForDel;
            statistics::InfoAboutContigsAlig infoAlig;
            int cnt[2] = {0, 0};
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void directSimpl(ContigGraph *graph, int v);

            void revSimpl(ContigGraph *graph, int v);

            int cntEdge(ContigGraph *graph, std::vector<int> edges);
        };
    }
}



#endif //SCAFFOLDER_RULECOORD_H

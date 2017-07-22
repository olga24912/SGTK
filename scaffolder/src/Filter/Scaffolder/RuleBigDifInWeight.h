#ifndef SCAFFOLDER_RULEBIGDIFINWEIGHT_H
#define SCAFFOLDER_RULEBIGDIFINWEIGHT_H

#include <ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleBigDifInWeight : public Rule {
        private:
            const int maxDif = 10;
        public:
            void simplifyGraph(ContigGraph *graph) override;

            void delSmallEdges(ContigGraph *graph, const std::vector<int> &edges) const;
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDIFINWEIGHT_H

#ifndef SCAFFOLDER_RULEBIGDIFINWEIGHT_H
#define SCAFFOLDER_RULEBIGDIFINWEIGHT_H

#include <ContigGraph/ContigGraph.h>
#include <set>
#include <Filter/Statistics/InfoAboutContigsAlig.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleBigDifInWeight : public Rule {
        private:
            std::set<int> edgeForDel;

        public:
            void simplifyGraph(ContigGraph *graph) override;

            void delSmallEdges(ContigGraph *graph, const std::vector<int> &edges);

            bool isGoodConnection(ContigGraph *graph, int e, const std::vector<int>& edges) const;

            bool isBadConnection(ContigGraph *graph, int e, const std::vector<int>& edges) const;
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDIFINWEIGHT_H

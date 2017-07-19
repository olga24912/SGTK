#ifndef SCAFFOLDER_RULEDELCYCLE_H
#define SCAFFOLDER_RULEDELCYCLE_H

#include <ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        class RuleDelCycle : public Rule {
        private:
            const int DEL_DIF = 2;

            std::vector<int> topsort;
            std::vector<int> color;
            std::vector<int> topSortPos;

            void topSortDfs(int v, ContigGraph *graph, std::vector<int> *used);

            void colorDfs(int v, int col, ContigGraph *graph);

            void topSort(ContigGraph *graph);

            int findCycle(ContigGraph *graph);
        public:
            void simplifyGraph(ContigGraph *graph) override;
        };
    }
}
#endif //SCAFFOLDER_RULEDELCYCLE_H

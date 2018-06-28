#ifndef SCAFFOLDER_RULEDELSMALLEDGES_H
#define SCAFFOLDER_RULEDELSMALLEDGES_H

#include <set>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDelSmallEdges : public Rule {
        private:
            const int startDel = 2;
            const int difDel = 3;

            std::set<int> delEdge;
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void delEdges(ContigGraph *graph, const std::vector<int>& edges);

        };
    }
}


#endif //SCAFFOLDER_RULEDELSMALLEDGES_H

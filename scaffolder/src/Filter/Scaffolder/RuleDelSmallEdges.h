#ifndef SCAFFOLDER_RULEDELSMALLEDGES_H
#define SCAFFOLDER_RULEDELSMALLEDGES_H

#include <set>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleDelSmallEdges : public Rule {
        private:
            const int startDel = 1;
            const int difDel = 2;

            std::set<int> delEdge;
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void delEdges(ContigGraph *graph, std::vector<int> edges);

            bool isAlone(ContigGraph *e, int i);

            bool sameCoord1(ContigGraph *graph, int e1, int e2);

            bool sameCoord2(ContigGraph *graph, int e1, int e2);
        };
    }
}


#endif //SCAFFOLDER_RULEDELSMALLEDGES_H

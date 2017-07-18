#ifndef SCAFFOLDER_RULEBIGDIFINWEIGHT_H
#define SCAFFOLDER_RULEBIGDIFINWEIGHT_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleBigDifInWeight : public Rule {
        private:
            const int maxDif = 5;
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void delSmallEdges(ContigGraph *filter, const std::vector<int> &edges) const;
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDIFINWEIGHT_H

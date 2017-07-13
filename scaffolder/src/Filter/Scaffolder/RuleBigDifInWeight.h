#ifndef SCAFFOLDER_RULEBIGDIFINWEIGHT_H
#define SCAFFOLDER_RULEBIGDIFINWEIGHT_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleBigDifInWeight : public Rule {
        private:
            const int maxDif = 3;
        public:
            void simplifyGraph(Filter *filter) override;

            void delSmallEdges(Filter *filter, const std::vector<int> &edges) const;
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDIFINWEIGHT_H

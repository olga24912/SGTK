#ifndef SCAFFOLDER_RULEBIGDEG_H
#define SCAFFOLDER_RULEBIGDEG_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleBigDeg : public Rule {
        private:
            const int BIG_DEG = 20;
            const int MAX_WEIGHT = 4;
        public:
            void simplifyGraph(Filter *filter) override;

            void ignoreEdges(Filter *wig, std::vector<int> vector);

            void ignoreEdgesR(Filter *wig, std::vector<int> vector);
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDEG_H

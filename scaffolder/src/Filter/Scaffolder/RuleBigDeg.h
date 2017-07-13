#ifndef SCAFFOLDER_RULEBIGDEG_H
#define SCAFFOLDER_RULEBIGDEG_H

#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleBigDeg : public Rule {
        private:
            const int BIG_DEG = 6;
            const int MAX_DIF = 2;
        public:
            void simplifyGraph(Filter *filter) override;

            void ignoreEdges(Filter *wig, std::vector<int> vector);

            void ignoreEdgesR(Filter *wig, std::vector<int> vector);
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDEG_H

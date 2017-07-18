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
            void simplifyGraph(ContigGraph *filter) override;

            void ignoreEdges(ContigGraph *wig, std::vector<int> vector);

            void ignoreEdgesR(ContigGraph *wig, std::vector<int> vector);
        };
    }
}


#endif //SCAFFOLDER_RULEBIGDEG_H

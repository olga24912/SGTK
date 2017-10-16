#ifndef SCAFFOLDER_RULEEXONBLOCKS_H
#define SCAFFOLDER_RULEEXONBLOCKS_H

#include <Filter/AligInfo/InfoAboutContigsAlig.h>
#include <Filter/Statistics/StrandStatistic.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        class RuleExonBlocks : public Rule {
        private:
            alig_info::InfoAboutContigsAlig infoAlig;
            filter::statistics::StrandStatistic* stat;
            int cnt[2] = {0, 0};
        public:
            void simplifyGraph(ContigGraph *filter) override;

            void checkEdge(ContigGraph *graph, int e);

            int checkFirstCoord(ContigGraph *graph, int e);

            void checkSecondCoord(ContigGraph *graph, int e);
        };
    }
}


#endif //SCAFFOLDER_RULEEXONBLOCKS_H

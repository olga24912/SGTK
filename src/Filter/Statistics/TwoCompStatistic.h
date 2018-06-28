#ifndef SCAFFOLDER_TWOCOMPSTATISTIC_H
#define SCAFFOLDER_TWOCOMPSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class TwoCompStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *filter, std::string coordFile,
                                    int libNum1, int val1,
                                    int libNum2, int val2);
        };
    }
}

#endif //SCAFFOLDER_TWOCOMPSTATISTIC_H

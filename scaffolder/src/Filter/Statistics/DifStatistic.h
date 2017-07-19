#ifndef SCAFFOLDER_DIFSTATISTIC_H
#define SCAFFOLDER_DIFSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class DifStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *filter, std::string coordFile, int libNum, int step, int mxVal);
        };
    }
}

#endif //SCAFFOLDER_DIFSTATISTIC_H

#ifndef SCAFFOLDER_WEIGHTDIFSTATISTIC_H
#define SCAFFOLDER_WEIGHTDIFSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class WeightDifStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *graph, std::string coordFile, int libNum, int stepW,
                                    int mxW, int stepDif, int mxDif);
        };
    }
}

#endif //SCAFFOLDER_WEIGHTDIFSTATISTIC_H

#ifndef SCAFFOLDER_WIEGHTSTATISTIC_H
#define SCAFFOLDER_WIEGHTSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class WeightStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *graph, std::string coordFile, int libNum, int step, int mxWeight);
        };
    }
}

#endif //SCAFFOLDER_WIEGHTSTATISTIC_H

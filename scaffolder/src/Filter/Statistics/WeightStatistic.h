#ifndef SCAFFOLDER_WIEGHTSTATISTIC_H
#define SCAFFOLDER_WIEGHTSTATISTIC_H

#include "Statistic.h"

namespace filter {
    namespace statistics {
        class WeightStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *graph, std::string coordFile, int libNum, int step, int mxWeight);
        };
    }
}

#endif //SCAFFOLDER_WIEGHTSTATISTIC_H

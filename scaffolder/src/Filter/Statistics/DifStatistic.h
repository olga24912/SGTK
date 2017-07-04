#ifndef SCAFFOLDER_DIFSTATISTIC_H
#define SCAFFOLDER_DIFSTATISTIC_H

#include "Statistic.h"

namespace filter {
    namespace statistics {
        class DifStatistic : public Statistic {
        public:
            void calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step, int mxVal);
        };
    }
}

#endif //SCAFFOLDER_DIFSTATISTIC_H

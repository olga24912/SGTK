#ifndef SCAFFOLDER_TWOLIBSTATISTIC_H
#define SCAFFOLDER_TWOLIBSTATISTIC_H

#include "Statistic.h"

namespace filter {
    namespace statistics {
        class TwoLibStatistic : public Statistic {
        public:
            void calculateStatistic(Filter *filter, std::string coordFile,
                                    int libNum1, int step1, int mxW1,
                                    int step1_2, int mwW1_2,
                                    int libNum2, int step2, int mxW2);
        };
    }
}
#endif //SCAFFOLDER_TWOLIBSTATISTIC_H

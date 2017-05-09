#ifndef SCAFFOLDER_TWOLIBSTATISTIC_H
#define SCAFFOLDER_TWOLIBSTATISTIC_H


#include "Statistic.h"

class TwoLibStatistic : public Statistic {
public:
    void calculateStatistic(Filter *filter, std::string coordFile,
                            int libNum1, int step1, int mxW1,
                            int libNum2, int step2, int mxW2);
};


#endif //SCAFFOLDER_TWOLIBSTATISTIC_H

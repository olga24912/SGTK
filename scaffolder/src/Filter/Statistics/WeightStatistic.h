#ifndef SCAFFOLDER_WIEGHTSTATISTIC_H
#define SCAFFOLDER_WIEGHTSTATISTIC_H


#include "Statistic.h"

class WeightStatistic : public Statistic {
public:
    void calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step, int mxWeight);
};


#endif //SCAFFOLDER_WIEGHTSTATISTIC_H

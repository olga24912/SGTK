#ifndef SCAFFOLDER_DIFSTATISTIC_H
#define SCAFFOLDER_DIFSTATISTIC_H

#include "Statistic.h"

class DifStatistic : public Statistic {
public:
    void calculateStatistic(Filter *filter, std::string coordFile, int libNum, int step, int mxVal);
};


#endif //SCAFFOLDER_DIFSTATISTIC_H

#ifndef SCAFFOLDER_TWOCOMPSTATISTIC_H
#define SCAFFOLDER_TWOCOMPSTATISTIC_H


#include "Statistic.h"

class TwoCompStatistic : public Statistic {
public:
    void calculateStatistic(Filter *filter, std::string coordFile,
                            int libNum1, int val1,
                            int libNum2, int val2);
};


#endif //SCAFFOLDER_TWOCOMPSTATISTIC_H

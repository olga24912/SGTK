#ifndef SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H


#include <Filter/Filters/Filter.h>
#include "Statistic.h"

class CntRightConnectionStatistic : public Statistic {
public:
    void calculateStatistic(Filter *filter, std::string coordFile, int libNum);

};


#endif //SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H

#ifndef SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H


#include <Filter/Filters/Filter.h>

class cntRightConnectionStatistic {
private:
    const int MAX_DIST = 1000000000;
public:
    void calculateStatistic(Filter *filter, std::string coordFile, int libNum);

};


#endif //SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H

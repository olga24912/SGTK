#ifndef SCAFFOLDER_WEIGHTHISTOGRAM_H
#define SCAFFOLDER_WEIGHTHISTOGRAM_H


#include <Filter/Filters/Filter.h>

class WeightHistogram {
public:
    static void histogram(Filter* filter, int lib, int step);
};


#endif //SCAFFOLDER_WEIGHTHISTOGRAM_H
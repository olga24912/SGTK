#ifndef SCAFFOLDER_STATISTIC_H
#define SCAFFOLDER_STATISTIC_H


#include "Filter/AligInfo/InfoAboutContigsAlig.h"

namespace filter {
    namespace statistics {
        class Statistic {
        protected:
            void printStatistic(char *str, int *cnt);
        };
    }
}


#endif //SCAFFOLDER_STATISTIC_H

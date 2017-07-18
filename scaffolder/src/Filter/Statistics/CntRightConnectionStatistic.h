#ifndef SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H


#include <Filter/Filters/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        class CntRightConnectionStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *filter, std::string coordFile, int libNum);

        };
    }
}

#endif //SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H

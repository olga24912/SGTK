#ifndef SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H

#include <ContigGraph/ContigGraph.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;

        class CntRightConnectionStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *filter, std::string coordFile, int libNum);

        };
    }
}

#endif //SCAFFOLDER_CNTRIGHTCONNECTIONSTATISTIC_H

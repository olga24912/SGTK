#ifndef SCAFFOLDER_STRANDSTATISTIC_H
#define SCAFFOLDER_STRANDSTATISTIC_H


#include <Filter/ContigGraph/ContigGraph.h>
#include <Filter/AligInfo/GeneAnnotationAns.h>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        using namespace contig_graph;

        class StrandStatistic : public Statistic {
        public:
            void calculateStatistic(ContigGraph *graph, std::string coordFile, std::string gffFile);
        };
    }
}



#endif //SCAFFOLDER_STRANDSTATISTIC_H

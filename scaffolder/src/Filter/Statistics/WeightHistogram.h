#ifndef SCAFFOLDER_WEIGHTHISTOGRAM_H
#define SCAFFOLDER_WEIGHTHISTOGRAM_H

#include <ContigGraph/ContigGraph.h>

namespace filter {
    namespace statistics {
        using namespace contig_graph;
        class WeightHistogram {
        public:
            static void histogram(ContigGraph *graph, int lib, int step);
        };
    }
}

#endif //SCAFFOLDER_WEIGHTHISTOGRAM_H

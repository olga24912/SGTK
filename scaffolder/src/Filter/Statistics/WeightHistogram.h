#ifndef SCAFFOLDER_WEIGHTHISTOGRAM_H
#define SCAFFOLDER_WEIGHTHISTOGRAM_H


#include <Filter/Filters/ContigGraph.h>

namespace filter {
    namespace statistics {
        class WeightHistogram {
        public:
            static void histogram(ContigGraph *graph, int lib, int step);
        };
    }
}

#endif //SCAFFOLDER_WEIGHTHISTOGRAM_H

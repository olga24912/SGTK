#ifndef SCAFFOLDER_SCAFFOLDERPIPELINE_H
#define SCAFFOLDER_SCAFFOLDERPIPELINE_H

#include <Filter/Filters/Filter.h>
#include <string>
#include "Scaffolds.h"

namespace filter {
    namespace scaffolder {
        class ScaffolderPipeline {
        private:

        public:
            void evaluate(Filter *graph, std::string contigFile, std::string out);

        private:
            DECL_LOGGER("ScaffolderPipeline");
        };
    }
}
#endif //SCAFFOLDER_SCAFFOLDERPIPELINE_H

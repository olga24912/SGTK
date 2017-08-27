#ifndef SCAFFOLDER_SCAFFOLDERPIPELINE_H
#define SCAFFOLDER_SCAFFOLDERPIPELINE_H

#include <string>
#include <ContigGraph/ContigGraph.h>
#include "Scaffolds.h"
#include <Filter/Statistics/InfoAboutContigsAlig.h>
#include <Filter/CommandParsers/State.h>

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        using namespace statistics;
        class ScaffolderPipeline {
        private:
            InfoAboutContigsAlig* alig;
        public:
            void evaluate(ContigGraph *graph, std::string contigFile,
                          std::string out, std::vector<commands::State::BamFiles> bamFiles);

            void setAlig(InfoAboutContigsAlig* alig) {
                this->alig = alig;
            }
        private:
            DECL_LOGGER("ScaffolderPipeline");
        };
    }
}
#endif //SCAFFOLDER_SCAFFOLDERPIPELINE_H

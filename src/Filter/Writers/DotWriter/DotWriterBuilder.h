#ifndef SCAFFOLDER_DOTWRITERBUILDER_H
#define SCAFFOLDER_DOTWRITERBUILDER_H


#include "DotWriter.h"

namespace filter {
    namespace writers {
        class DotWriterBuilder {
        protected:
            ContigGraph *graph;
            FileValidator *validator;
            std::string coordFile = "";
            int maxVert = 20;
            int maxEdge = 40;
        public:
            virtual DotWriter *build() {
                DEBUG("build simple dot writer");
                return new DotWriter(graph, validator, maxVert, maxEdge, coordFile);
            }

            void setFilter(ContigGraph *filter) {
                this->graph = filter;
            }

            void setValidator(FileValidator *validator) {
                this->validator = validator;
            }

            void setMaxVert(int maxVert) {
                DotWriterBuilder::maxVert = maxVert;
            }

            void setMaxEdge(int maxEdge) {
                DotWriterBuilder::maxEdge = maxEdge;
            }

            void setCoordFile(std::string coordFile) {
                DotWriterBuilder::coordFile = coordFile;
            }
        protected:
            DECL_LOGGER("DotWriterBuilder");
        };
    }
}

#endif //SCAFFOLDER_DOTWRITERBUILDER_H

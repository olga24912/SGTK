#ifndef SCAFFOLDER_CHRDOTWRITERBUILDER_H
#define SCAFFOLDER_CHRDOTWRITERBUILDER_H

#include "DotWriterBuilder.h"
#include "ChrDotWriter.h"

namespace filter {
    namespace writers {
        class ChrDotWriterBuilder : public DotWriterBuilder {
        public:
            virtual DotWriter *build() override {
                DEBUG("build chr DW");
                return new ChrDotWriter(graph, validator, coordFile);
            }
        };
    }
}

#endif //SCAFFOLDER_CHRDOTWRITERBUILDER_H

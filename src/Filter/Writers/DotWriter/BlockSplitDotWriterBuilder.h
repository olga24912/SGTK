#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

#include "DotWriterBuilder.h"
#include "BlockSplitDotWriter.h"

namespace filter {
    namespace writers {
        class BlockSplitDotWriterBuilder : public DotWriterBuilder {
        public:
            virtual DotWriter *build() override {
                DEBUG("build block DW");
                return new BlockSplitDotWriter(graph, validator, maxVert, maxEdge);
            }
        };
    }
}

#endif //SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H


#include "DotWriterBuilder.h"
#include "BlockSplitDotWriter.h"

class BlockSplitDotWriterBuilder : public DotWriterBuilder {
public:
    DotWriter build() override {
        return BlockSplitDotWriter(filter, validator, maxVert, maxEdge);
    }
};

#endif //SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

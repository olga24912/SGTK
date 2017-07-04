#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H


#include "DotWriterBuilder.h"
#include "BlockSplitDotWriter.h"

class BlockSplitDotWriterBuilder : public DotWriterBuilder {
public:
    virtual DotWriter* build() override {
        DEBUG("build block DW");
        return new BlockSplitDotWriter(filter, validator, maxVert, maxEdge);
    }
};

#endif //SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

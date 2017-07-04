#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H


#include "DotWriterBuilder.h"
#include "BlockSplitDotWriter.h"

class BlockSplitDotWriterBuilder : public DotWriterBuilder {
public:
    virtual DotWriter* build() override {

        std::cerr << "build block DW" << std::endl;
        return new BlockSplitDotWriter(filter, validator, maxVert, maxEdge);
    }
};

#endif //SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

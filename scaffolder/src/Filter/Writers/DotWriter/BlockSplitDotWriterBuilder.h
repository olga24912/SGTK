#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H


#include "DotWriterBuilder.h"

class BlockSplitDotWriterBuilder : public DotWriterBuilder {
public:
    DotWriter build() override {
        return DotWriterBuilder::build();
    }
};

#endif //SCAFFOLDER_BLOCKSPLITDOTWRITERBUILDER_H

//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H
#define SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

#include "GraphBuilder.h"

class RNASplitReadGraphBuilder : public GraphBuilder {
private:
    string refFileName;
    string rnaReadsFileName;

public:
    void evaluate();

    void setRefFileName(string refFileName);

    void setRnaReadFileName(string rnaReadsFileName);
};


#endif //SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

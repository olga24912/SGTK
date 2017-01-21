//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

class RNAPairReadGraphBuilder : public PairReadGraphBuilder {
private:
    virtual string getLibColor();
public:
    RNAPairReadGraphBuilder(){}
};


#endif //SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

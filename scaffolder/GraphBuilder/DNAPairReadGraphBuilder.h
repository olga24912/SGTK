//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

class DNAPairReadGraphBuilder: public PairReadGraphBuilder {
private:
    int distBetweenPairReads;

public:
    void setDistBetweenPairReads(int distBetweenPairReads);

};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

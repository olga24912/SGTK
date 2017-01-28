#ifndef SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

//generate connection between contigs by RNA pair reads
class RNAPairReadGraphBuilder : public PairReadGraphBuilder {
private:
    virtual std::string getLibColor();
public:
    RNAPairReadGraphBuilder(){}
};


#endif //SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

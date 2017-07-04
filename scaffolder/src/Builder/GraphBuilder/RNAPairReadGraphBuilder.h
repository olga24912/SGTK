#ifndef SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

//generate connection between contigs by RNA pair reads
class RNAPairReadGraphBuilder : public PairReadGraphBuilder {
private:
    virtual std::string getLibColor();

protected:
    ContigGraph::Lib::Type getLibType() override;

public:
    RNAPairReadGraphBuilder(){}
private:
    DECL_LOGGER("RNAPairReadGraphBuilder");
};


#endif //SCAFFOLDER_RNAPAIRREADGRAPHBUILDER_H

//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

class DNAPairReadGraphBuilder: public PairReadGraphBuilder {
private:
    int distBetweenPairReads;
    map<string, int> read1DistToEnd;
    map<string, int> read2DistToEnd;

public:
    void setDistBetweenPairReads(int distBetweenPairReads);

    void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);

    int readDist(BamAlignmentRecord read);
};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

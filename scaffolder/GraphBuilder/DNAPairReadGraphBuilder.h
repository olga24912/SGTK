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

    void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);

    int readDist(BamAlignmentRecord read);
public:
    void setDistBetweenPairReads(int distBetweenPairReads);

};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"

class DNAPairReadGraphBuilder: public PairReadGraphBuilder {
private:
    static constexpr double DEFAULT_DEF = 0.179;
    static const int DEFAULT_MAX_COUNT_EDGE = 3;

    int distBetweenPairReads;
    map<string, int> read1DistToEnd;
    map<string, int> read2DistToEnd;

    void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);

    void incEdgeWeight(string readName, int target);

    void filterEdge();

    int readDist(BamAlignmentRecord read);
public:
    void setDistBetweenPairReads(int distBetweenPairReads);

    int countEdgesBeforeBreak(int v, vector<int> edges);
};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

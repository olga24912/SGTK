//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

#include "PairReadGraphBuilder.h"
/*
 * find conection between contigs by DNA pair reads
 */
class DNAPairReadGraphBuilder: public PairReadGraphBuilder {
private:
    static constexpr double DEFAULT_DEF = 0.179;
    static const int DEFAULT_MAX_COUNT_EDGE = 3;

    int distBetweenPairReads;
    map<string, int> read1DistToEnd;
    map<string, int> read2DistToEnd;

    void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);

    void incEdgeWeight(BamAlignmentRecord read1, BamAlignmentRecord read2);

    void filterEdge();

    int readDist(BamAlignmentRecord read);

    int countEdgesBeforeBreak(int v, vector<int> edges);
public:
    /*
     * set max distance between DNA pair read,
     * if  dist will be more this conatcion will be ignore/
     */
    void setDistBetweenPairReads(int distBetweenPairReads);
};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

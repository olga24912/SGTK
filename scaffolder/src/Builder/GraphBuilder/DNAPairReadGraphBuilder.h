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
    int distBetweenPairReads;
    map<string, int> read1DistToEnd;
    map<string, int> read2DistToEnd;

    void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);

    void incEdgeWeight(BamAlignmentRecord read1, BamAlignmentRecord read2);

    int readDist(BamAlignmentRecord read);

    virtual string getLibColor();
public:
    /*
     * set max distance between DNA pair read,
     * if  dist will be more this conatcion will be ignore/
     */
    void setDistBetweenPairReads(int distBetweenPairReads);
};


#endif //SCAFFOLDER_DNAPAIRREADGRAPHBUILDER_H

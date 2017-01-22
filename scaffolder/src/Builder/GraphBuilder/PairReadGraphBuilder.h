//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_PAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_PAIRREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include "SamFileWriter/SamFileWriter.h"
#include "SamFileWriter/SamFileWriteEdge.h"
#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>


using namespace std;
using namespace seqan;

//genrate conection betwen contigs by pair reads
class PairReadGraphBuilder: public GraphBuilder {
protected:
    bool oneSideRead = false;

    string fileName1;
    string fileName2;

    BamFileIn bamFile1;
    BamFileIn bamFile2;

    unordered_map<string, BamAlignmentRecord> read1ByName;
    unordered_map<string, BamAlignmentRecord> read2ByName;

    pair<string, int> processOneFirstRead(BamAlignmentRecord read);
    pair<string, int> processOneSecondRead(BamAlignmentRecord read);
    virtual void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    virtual void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);
    int get2Target(const BamAlignmentRecord &read) const;
    int get1Target(const BamAlignmentRecord &read) const;
    void readHeaderInit();
    void addInfoAboutCover(int target, const BamAlignmentRecord &read);
    bool isUniqueMapRead(BamAlignmentRecord read);

    virtual void incEdgeWeight(BamAlignmentRecord read1, BamAlignmentRecord read2);
    virtual void filterEdge();

    int pairTarget(int id);

    void handleReads();
public:
    //set sam file for first pair read alignment
    void setFileName2(const string &fileName2);

    //set sam file for second pair read alignment
    void setFileName1(const string &fileName1);

    //set info about orintation  pair read
    void setOneSideReadFlag(bool flag);

    virtual void evaluate();
};


#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

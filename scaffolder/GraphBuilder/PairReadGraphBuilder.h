//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_PAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_PAIRREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>


using namespace std;
using namespace seqan;

class PairReadGraphBuilder: public GraphBuilder {
protected:
    bool oneSideRead = false;

    string fileName1;
    string fileName2;

    BamFileIn bamFile1;
    BamFileIn bamFile2;

    unordered_map<string, int> read1Target;
    unordered_map<string, int> read2Target;

    pair<string, int> processOneFirstRead(BamAlignmentRecord read);
    pair<string, int> processOneSecondRead(BamAlignmentRecord read);
    virtual void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
    virtual void addInfoAbout2Read(string readName, int target, BamAlignmentRecord read);
    void readHeaderInit();
    void addInfoAboutCover(int target, const BamAlignmentRecord &read);

    virtual void incEdgeWeight(string readName, int target);
    virtual void filterEdge();

    int pairTarget(int id);
public:
    void setFileName2(const string &fileName2);
    void setFileName1(const string &fileName1);

    void setOneSideReadFlag(bool flag);

    virtual void evaluate();

    void handleReads();
};


#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

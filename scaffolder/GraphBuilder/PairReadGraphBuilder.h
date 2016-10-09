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
    string fileName1;
    string fileName2;

    BamFileIn bamFile;

    map<string, int> read1Target;

    void firstReads();
    void secondReads();

    void processOneFirstRead(BamAlignmentRecord read);
    virtual void addInfoAboutRead(string readName, int target, BamAlignmentRecord read);
public:
    void setFileName2(const string &fileName2);
    void setFileName1(const string &fileName1);

    virtual void evaluate();

    void readHeaderInit();
};


#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

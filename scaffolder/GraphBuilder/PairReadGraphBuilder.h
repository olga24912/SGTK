//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_PAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_PAIRREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include <bits/stdc++.h>

using namespace std;
using namespace seqan;

class PairReadGraphBuilder: public GraphBuilder {
protected:
    string fileName1;
    string fileName2;

    BamFileIn bamFile;

    void firstReads();
    void secondReads();
public:
    void setFileName2(const string &fileName2);
    void setFileName1(const string &fileName1);

    virtual void evaluate();

    void readHeaderInit();
};


#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

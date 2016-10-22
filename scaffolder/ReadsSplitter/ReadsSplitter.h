//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_READSSPLITER_H
#define SCAFFOLDER_READSSPLITER_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include "../Tools/FastaToolsOut.h"


using namespace seqan;
using namespace std;

class ReadsSplitter {
public:
    unordered_map<string, string> getFullReads() const {
        return reads;
    }
protected:
    unordered_map<string, string> reads;

    void splitRead(string readName, string seq, int len, FastaToolsOut& out1, FastaToolsOut& out2);
public:
    virtual void splitReads(string rnaFileName, string resFileName1, string resFileName2) = 0;
};


#endif //SCAFFOLDER_READSSPLITER_H

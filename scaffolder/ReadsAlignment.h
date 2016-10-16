//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_READSALIGNMENT_H
#define SCAFFOLDER_READSALIGNMENT_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>
#include "Tools/SeqanUtils.h"



using namespace seqan;
using namespace std;

class ReadsAlignment {
public:
    void alignment(string refFileName, unordered_map<string, string> fullReads,
                   string alignfileName, string resFile1, string resFile2);

private:
    const int inf = 1e9;

    unordered_map<string, string> getContigs(string fileName);

    int alignFull(string ref, string read, int pos);
};


#endif //SCAFFOLDER_READSALIGNMENT_H

//
// Created by olga on 31.10.16.
//

#ifndef SCAFFOLDER_SAMFILEWRITEEDGE_H
#define SCAFFOLDER_SAMFILEWRITEEDGE_H

#include "bits/stdc++.h"
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

using namespace std;
using namespace seqan;

class SamFileWriteEdge {
private:
    string dir;
    string getName(int edgeID);
    BamFileIn* fileIn;

public:
    SamFileWriteEdge(string dir);
    SamFileWriteEdge(){};
    void setFileIn(BamFileIn* in) {
        fileIn = in;
    }
    virtual void writeEdge(int edgeID, BamAlignmentRecord read1, BamAlignmentRecord read2);
};


#endif //SCAFFOLDER_SAMFILEWRITEEDGE_H

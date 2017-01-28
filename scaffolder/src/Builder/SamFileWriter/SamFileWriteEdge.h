//
// Created by olga on 31.10.16.
//

#ifndef SCAFFOLDER_SAMFILEWRITEEDGE_H
#define SCAFFOLDER_SAMFILEWRITEEDGE_H

#include "bits/stdc++.h"
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

class SamFileWriteEdge {
private:
    std::string dir;
    std::string getName(int edgeID);
    seqan::BamFileIn* fileIn;

public:
    SamFileWriteEdge(std::string dir);
    SamFileWriteEdge(){};
    void setFileIn(seqan::BamFileIn* in) {
        fileIn = in;
    }
    virtual void writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);
};


#endif //SCAFFOLDER_SAMFILEWRITEEDGE_H

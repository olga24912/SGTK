//
// Created by olga on 31.10.16.
//

#ifndef SCAFFOLDER_SAMFILEWRITEEDGE_H
#define SCAFFOLDER_SAMFILEWRITEEDGE_H

#include "SamFileWriter.h"
#include "bits/stdc++.h"

using namespace std;

class SamFileWriteEdge : public SamFileWriter {
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

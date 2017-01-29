#ifndef SCAFFOLDER_SAMFILEWRITEEDGE_H
#define SCAFFOLDER_SAMFILEWRITEEDGE_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

//write info about reads for current edge
class SamFileWriteEdge {
private:
    std::string dir;
    std::string getName(int edgeID);
    seqan::BamFileIn* fileIn;

public:
    //dir where fails will be create
    SamFileWriteEdge(std::string dir);
    SamFileWriteEdge(){};
    //need for context for create BamFileOut
    void setFileIn(seqan::BamFileIn* in) {
        fileIn = in;
    }
    //write record for edgeId with read1 and read2
    virtual void writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);
};


#endif //SCAFFOLDER_SAMFILEWRITEEDGE_H

#ifndef SCAFFOLDER_SAMFILEWRITEREDGEFICTIVE_H
#define SCAFFOLDER_SAMFILEWRITEREDGEFICTIVE_H

#include "SamFileWriteEdge.h"

class SamFileWriterEdgeFictive : public SamFileWriteEdge {
public:
    SamFileWriterEdgeFictive(std::string dir);

    SamFileWriterEdgeFictive();

    void setFileIn(seqan::BamFileIn *in) override;

    void writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) override;
};


#endif //SCAFFOLDER_SAMFILEWRITEREDGEFICTIVE_H

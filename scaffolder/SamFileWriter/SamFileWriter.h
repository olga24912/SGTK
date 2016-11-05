//
// Created by olga on 31.10.16.
//

#ifndef SCAFFOLDER_SAMFILEWRITER_H
#define SCAFFOLDER_SAMFILEWRITER_H

#include <seqan/bam_io.h>

using namespace seqan;

class SamFileWriter {
public:
    virtual void writeEdge(int edgeID, BamAlignmentRecord read1, BamAlignmentRecord read2){};
};


#endif //SCAFFOLDER_SAMFILEWRITER_H

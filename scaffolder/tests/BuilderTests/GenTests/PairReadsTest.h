#ifndef SCAFFOLDER_PAIRREADSTEST_H
#define SCAFFOLDER_PAIRREADSTEST_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "Builder/Tools/SeqanUtils.h"

class PairReadsTest {
protected:
    int pairDist = 200;
    int readLen = 100;
    int contigLen = 1000;

    std::string dnaSeq;

    std::string read1FileName, read2FileName;
    std::string outRefFileNema;
    std::string refFileName;

    void genPairRead(std::string ref);
    void genContig();
    std::string rev_compl(std::string s);

public:
    virtual void genTest(std::string refFileName, std::string read1FileName,
                 std::string read2FileName, std::string outRefFileName) = 0;

};


#endif //SCAFFOLDER_PAIRREADSTEST_H

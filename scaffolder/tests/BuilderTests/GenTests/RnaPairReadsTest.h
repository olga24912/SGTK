#ifndef SCAFFOLDER_RNAPAIRREADSTEST_H
#define SCAFFOLDER_RNAPAIRREADSTEST_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "Builder/Tools/SeqanUtils.h"

class RnaPairReadsTest {
private:
    int exonLen = 1500;
    int intrLen = 500;
    int pairDist = 200;
    int readLen = 100;
    int contigLen = 1000;

    std::string dnaSeq;
    std::string rnaSeq;

    std::string read1FileName, read2FileName;
    std::string outRefFileNema;
    std::string refFileName;

    void genPairRead();
    void genRnaSeq();
    void genContig();
    std::string rev_compl(std::string s);

public:
    void genTest(std::string refFileName, std::string read1FileName,
                 std::string read2FileName, std::string outRefFileName);

};


#endif //SCAFFOLDER_RNAPAIRREADSTEST_H

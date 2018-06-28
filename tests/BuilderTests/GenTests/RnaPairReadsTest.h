#ifndef SCAFFOLDER_RNAPAIRREADSTEST_H
#define SCAFFOLDER_RNAPAIRREADSTEST_H

#include "PairReadsTest.h"

class RnaPairReadsTest : public PairReadsTest {
private:
    int exonLen = 1500;
    int intrLen = 500;

    std::string rnaSeq;

    void genRnaSeq();

public:
    void genTest(std::string refFileName, std::string read1FileName,
                 std::string read2FileName, std::string outRefFileName);

};


#endif //SCAFFOLDER_RNAPAIRREADSTEST_H

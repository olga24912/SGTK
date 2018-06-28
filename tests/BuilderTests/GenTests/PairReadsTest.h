#ifndef SCAFFOLDER_PAIRREADSTEST_H
#define SCAFFOLDER_PAIRREADSTEST_H

#include "GenTest.h"

class PairReadsTest : public GenTest {
protected:
    int pairDist = 200;
    int readLen = 100;

    std::string read1FileName, read2FileName;

    void genPairRead(std::string ref);

public:
    virtual void genTest(std::string refFileName, std::string read1FileName,
                 std::string read2FileName, std::string outRefFileName) = 0;

};


#endif //SCAFFOLDER_PAIRREADSTEST_H

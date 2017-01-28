#ifndef SCAFFOLDER_DNAPAIRREADSTEST_H
#define SCAFFOLDER_DNAPAIRREADSTEST_H

#include "PairReadsTest.h"

class DnaPairReadsTest : public PairReadsTest {
public:
    void genTest(std::string refFileName, std::string read1FileName, std::string read2FileName,
                 std::string outRefFileName) override;
};


#endif //SCAFFOLDER_DNAPAIRREADSTEST_H

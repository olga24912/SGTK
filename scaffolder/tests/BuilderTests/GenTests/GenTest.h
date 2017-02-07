#ifndef SCAFFOLDER_REFGENTEST_H
#define SCAFFOLDER_REFGENTEST_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "ContigGraph/SeqanUtils.h"

class GenTest {
protected:
    int contigLen = 1000;

    std::string outRefFileNema;
    std::string refFileName;

    std::string dnaSeq;

    void genContig();
    std::string rev_compl(std::string s);

public:
    void genTest(std::string refFileName, std::string outRefFileName);
};


#endif //SCAFFOLDER_REFGENTEST_H

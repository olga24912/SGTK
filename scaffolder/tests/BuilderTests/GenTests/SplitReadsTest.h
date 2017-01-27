#ifndef SCAFFOLDER_SPLITREADSTEST_H
#define SCAFFOLDER_SPLITREADSTEST_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "Builder/Tools/SeqanUtils.h"

/*
 * Class for generate small test for slpit read
 * Take file with one ref genome if fasta format. Split it on two contigs
 * generate reads without mistake from some merge part of two contigs;
 */
class SplitReadsTest {
private:
    std::string smallRef(std::string ref, int exonLen, int intronLen, std::vector<int> spos);
public:
    /*
     * function  for generate test
     * fi - input file in fasta format with ref gen.
     * fo - output file for ref contigs.
     * rf - output file for generate reads.
     *
     * readLen - len of generate reads.
     */
    void genTest(std::string fi, std::string fo, std::string rf, int readLen, int contigCnt=2);
};


#endif //SCAFFOLDER_SPLITREADSTEST_H

//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_SPLITREADSTEST_H
#define SCAFFOLDER_SPLITREADSTEST_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "../Tools/SeqanUtils.h"

using namespace seqan;
using namespace std;

/*
 * Class for generate small test for slpit read
 * Take file with one ref genome if fasta format. Split it on two contigs
 * generate reads without mistake from some merge part of two contigs;
 */
class SplitReadsTest {
private:
    string smallRef(string ref, int exonLen, int intronLen, int spos);
public:
    /*
     * function  for generate test
     * fi - input file in fasta format with ref gen.
     * fo - output file for two ref contigs.
     * rf - output file for generate reads.
     *
     * readLen - len of generate reads.
     */
    void genTest(string fi, string fo, string rf, int readLen);
};


#endif //SCAFFOLDER_SPLITREADSTEST_H

//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_READSSPLITER_H
#define SCAFFOLDER_READSSPLITER_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>


using namespace seqan;
using namespace std;

/*
 * take reads from fasta file and split it on two equal parts.
 */
class ReadsSplitter {
public:
    unordered_map<string, string> getFullNotAlignmentReads() const {
        return reads;
    }
private:
    unordered_map<string, string> reads;
public:
    /*
     * take read from rnaUnmappedReadsFileName
     * first part of split read will be in resFileName1
     * second - resFileName2
     */
    void splitNotAlignmentReads(string rnaUnmappedReadsFileName,
                                string resFileName1, string resFileName2);
};


#endif //SCAFFOLDER_READSSPLITER_H

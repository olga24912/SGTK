//
// Created by olga on 05.11.16.
//

#ifndef SCAFFOLDER_CONTIGMERGER_H
#define SCAFFOLDER_CONTIGMERGER_H

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include "../Tools/SeqanUtils.h"
#include"bits/stdc++.h"

using namespace std;
using namespace seqan;

class ContigMerger {
private:
    int countN = 30;

    string contig1Name;
    string contig2Name;

    string contig1Val;
    string contig2Val;

    string contigName = "contig";
    string contigVal;

    string findContig(string fileIn, string name);
    void mergeContigs();
    void writeContig(string fileName);
    void writeHeader(BamFileOut& out);
public:
    void evaluate(string contigsINFileName, string samReads1FileName,
                  string samReads2FileName, string contigOUTFileName,
                  string samOutFileName,
                  string contig1Name, string contig2Name);
};


#endif //SCAFFOLDER_CONTIGMERGER_H

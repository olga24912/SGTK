//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_SEQANUTILS_H
#define SCAFFOLDER_SEQANUTILS_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

//handling info from seqans classes
class SeqanUtils {
public:
    static string cutReadName(BamAlignmentRecord read); //return reads name without "/1", "/2" end.
    static string dna5ToString(Dna5* seq, int len);
    static void writeRec(SeqFileOut& out, string name, string seq);
};


#endif //SCAFFOLDER_SEQANUTILS_H

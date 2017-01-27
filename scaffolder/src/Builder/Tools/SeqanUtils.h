//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_SEQANUTILS_H
#define SCAFFOLDER_SEQANUTILS_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

//handling info from seqans classes
class SeqanUtils {
public:
    static std::string cutReadName(seqan::BamAlignmentRecord read); //return reads name without "/1", "/2" end.
    static std::string dna5ToString(seqan::Dna5* seq, int len);
    static void writeRec(seqan::SeqFileOut& out, std::string name, std::string seq);
};


#endif //SCAFFOLDER_SEQANUTILS_H

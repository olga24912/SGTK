//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_SEQANUTILS_H
#define SCAFFOLDER_SEQANUTILS_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>


using namespace seqan;
using namespace std;

/*
 * handling info from seqans classes
 */
class SeqanUtils {
public:
    /*
     * return reads name without "/1", "/2" end.
     */
    static string cutReadName(BamAlignmentRecord read);
};


#endif //SCAFFOLDER_SEQANUTILS_H

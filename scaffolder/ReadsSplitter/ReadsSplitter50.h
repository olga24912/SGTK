//
// Created by olga on 22.10.16.
//

#ifndef SCAFFOLDER_READSSPLITTER50_H
#define SCAFFOLDER_READSSPLITTER50_H


#include "ReadsSplitter.h"
#include "../Tools/SeqanUtils.h"
#include <seqan/seq_io.h>

using namespace std;
using namespace seqan;

class ReadsSplitter50 : public ReadsSplitter {
public:
    /*
     * take read from rnaFileName
     * first part of split read will be in resFileName1
     * second - resFileName2
     */
    virtual void splitReads(string rnaFileName,
                            string resFileName1, string resFileName2);

};


#endif //SCAFFOLDER_READSSPLITTER50_H

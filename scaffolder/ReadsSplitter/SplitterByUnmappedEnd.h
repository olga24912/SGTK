//
// Created by olga on 22.10.16.
//

#ifndef SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H
#define SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H


#include "ReadsSplitter.h"
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

using namespace std;
using namespace seqan;

class SplitterByUnmappedEnd : public ReadsSplitter {
private:
    const int MIN_READ_LEN = 10;
public:
    virtual void splitReads(string rnaFileName,
                            string resFileName1, string resFileName2);
};


#endif //SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H

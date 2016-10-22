//
// Created by olga on 22.10.16.
//

#include "ReadsSplitter50.h"

void ReadsSplitter50::splitReads(string rnaUnmappedReadsFileName, string resFileName1, string resFileName2) {
    cerr << "start split reads" << endl;
    FastaToolsIn ftin;
    ftin.parse(rnaUnmappedReadsFileName);

    FastaToolsOut ftout1;
    ftout1.putFileName(resFileName1);

    FastaToolsOut ftout2;
    ftout2.putFileName(resFileName2);

    cerr << "start rewrite reads" << endl;

    while (ftin.next()) {
        string readName = ftin.currentName();
        string seq = ftin.currentRef();
        reads[readName] = seq;

        int len = (int) seq.size() / 2;
        splitRead(readName, seq, len, ftout1, ftout2);
    }

    ftin.close();
    ftout1.close();
    ftout2.close();
}


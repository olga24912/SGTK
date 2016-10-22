//
// Created by olga on 16.10.16.
//

#include "ReadsSplitter.h"
#include "../Tools/FastaToolsIn.h"
#include "../Tools/FastaToolsOut.h"


void ReadsSplitter::splitRead(string readName, string seq, int len, FastaToolsOut& out1, FastaToolsOut& out2) {
    string readName1 = readName;
    readName1 += "/1";
    string readName2 = readName;
    readName2 += "/2";

    out1.write(readName1, seq.substr(0, len));
    out2.write(readName2, seq.substr(len, seq.size() - len));
}

//
// Created by olga on 10.10.16.
//

#ifndef SCAFFOLDER_FASTATOOLS_H
#define SCAFFOLDER_FASTATOOLS_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>

using namespace std;
using namespace seqan;

class FastaToolsIn {
private:
    string fileName;
    ifstream fin;

    string curHeader;
    string curGenRef;
public:
    void parse(string fn);
    bool next();
    string currentName();
    string currentRef();
};


#endif //SCAFFOLDER_FASTATOOLS_H

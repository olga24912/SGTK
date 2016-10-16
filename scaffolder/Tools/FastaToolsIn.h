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
    int pos;
    ifstream fin;
public:
    void parse(string fn);
    void next();
    void currentName();
    void currentRef();
};


#endif //SCAFFOLDER_FASTATOOLS_H

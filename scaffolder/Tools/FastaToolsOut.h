//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_FASTATOOLSOUT_H
#define SCAFFOLDER_FASTATOOLSOUT_H

#include <bits/stdc++.h>

using namespace std;

class FastaToolsOut {
private:
    string fileName;

    ofstream fout;
public:
    void putFileName(string fileName);
    void write(string name, string ref);
    void close();
};


#endif //SCAFFOLDER_FASTATOOLSOUT_H

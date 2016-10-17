//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_FASTATOOLSOUT_H
#define SCAFFOLDER_FASTATOOLSOUT_H

#include <bits/stdc++.h>

using namespace std;
/*
 * class for print info to fasta format.
 */
class FastaToolsOut {
private:
    string fileName;

    ofstream fout;
public:
    /*
     * put file, that we going to write.
     */
    void putFileName(string fileName);

    /*
     * write info about ref with this name and this value
     */
    void write(string name, string ref);

    /*
     * need to call after work. Close output stream.
     */
    void close();
};


#endif //SCAFFOLDER_FASTATOOLSOUT_H

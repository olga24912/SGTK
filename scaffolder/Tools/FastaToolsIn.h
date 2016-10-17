//
// Created by olga on 10.10.16.
//

#ifndef SCAFFOLDER_FASTATOOLS_H
#define SCAFFOLDER_FASTATOOLS_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>

using namespace std;
using namespace seqan;

/*
 * class for reading from fasta file
 */
class FastaToolsIn {
private:
    string fileName;
    ifstream fin;

    string curHeader;
    string curGenRef;
    string nextHeader;
public:
    /*
     * set the file name in fasta format witch from we gone to read.
     */
    void parse(string fn);

    /*
     * read next pattern include name and gene
     */
    bool next();

    /*
     * get the name of read/contig
     */
    string currentName();

    /*
     * get gens value
     */
    string currentRef();

    /*
     * need to call after work. Close input stream.
     */
    void close();
};


#endif //SCAFFOLDER_FASTATOOLS_H

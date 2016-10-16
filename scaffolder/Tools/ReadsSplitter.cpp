//
// Created by olga on 16.10.16.
//

#include "ReadsSplitter.h"

void ReadsSplitter::findAndSplitNotAlignmentReads(string rnaReadsFileName, string alignmentFileName, string resFileName) {
    cerr << "start split reads" << endl;
    unordered_set<string> alignmentReads = findAlignmentReads(alignmentFileName);

    ofstream fout(resFileName);
    ifstream fin(rnaReadsFileName);

    cerr << "start rewrite reads" << endl;
    string cur;
    getline(fin, cur);
    while (true) {
        string readName = getName(cur);

        if (alignmentReads.count(readName) == 0) {
            string seq;
            getline(fin, seq);

            int len = (int)seq.size()/2;
            string readName1 = readName;
            readName1 += "_1";
            string readName2 = readName;
            readName2 += "_2";

            string next;
            if (getline(fin, next)) {
                if (next[0] == '+') {
                    string nn;
                    getline(fin, nn);

                    fout << cur[0] << readName1 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                    fout << seq.substr(0, len) << "\n";
                    fout << "+\n";
                    fout << nn.substr(0, len) << "\n";

                    fout << cur[0] << readName2 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                    fout << seq.substr(len, seq.size() - len) << "\n";
                    fout << "+\n";
                    fout << nn.substr(len, seq.size() - len) << "\n";

                    if (!getline(fin, cur)) {
                        break;
                    }

                } else {
                    fout << cur[0] << readName1 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                    fout << seq.substr(0, len) << "\n";

                    fout << cur[0] << readName2 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                    fout << seq.substr(len, seq.size() - len) << "\n";
                    cur = next;
                }
            } else {
                fout << cur[0] << readName1 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                fout << seq.substr(0, len) << "\n";

                fout << cur[0] << readName2 << cur.substr(readName.size() + 1, cur.size() - readName.size() - 1) << "\n";
                fout << seq.substr(len, seq.size() - len) << "\n";
                break;
            }
        } else {
            if (!getline(fin, cur)) {
                break;
            }
        }
    }

    fout.close();
    fin.close();
}

unordered_set<string> ReadsSplitter::findAlignmentReads(string fileName) {
    cerr << "start read sam file" << endl;
    unordered_set<string> usedReads;

    BamFileIn bamFile;
    open(bamFile, fileName.c_str());

    BamHeader samHeader;
    readHeader(samHeader, bamFile);

    BamAlignmentRecord read;
    while (!atEnd(bamFile)) {
        readRecord(read, bamFile);
        string name = string(toCString(read.qName));
        cerr << name << "\n";
        if (read.rID != -1) {
            usedReads.insert(name);
        }
    }
    close(bamFile);
    return usedReads;
}

string ReadsSplitter::getName(string headerLine) {
    string name;
    for (int i = 1; i < (int)headerLine.size() && headerLine[i] != ' '; ++i) {
        name += headerLine[i];
    }

    return name;
}
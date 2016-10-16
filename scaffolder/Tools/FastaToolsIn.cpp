//
// Created by olga on 10.10.16.
//

#include "FastaToolsIn.h"


void FastaToolsIn::parse(string fn) {
    fileName = fn;
    fin.open(fn);
}

bool FastaToolsIn::next() {
    if (getline(fin, curHeader)) {
        getline(fin, curGenRef);
        return true;
    }
    return false;
}

string FastaToolsIn::currentName() {
    string name;
    for (int i = 1; i < (int)curHeader.size() && curHeader[i] != ' '; ++i) {
        name += curHeader[i];
    }

    return name;
}


string FastaToolsIn::currentRef() {
    return curGenRef;
}

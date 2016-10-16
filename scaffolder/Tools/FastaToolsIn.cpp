//
// Created by olga on 10.10.16.
//

#include "FastaToolsIn.h"


void FastaToolsIn::parse(string fn) {
    fileName = fn;
    fin.open(fn);

    getline(fin, nextHeader);
}

bool FastaToolsIn::next() {
    curHeader = nextHeader;
    if (curHeader == "") return false;
    nextHeader = "";
    curGenRef = "";
    string s;
    while (getline(fin, s) && s[0] != '>') {
        curGenRef += s;
    }
    if (s.size() > 0 && s[0] == '>') {
        nextHeader = s;
    }
    return true;
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

void FastaToolsIn::close() {
    fin.close();
}

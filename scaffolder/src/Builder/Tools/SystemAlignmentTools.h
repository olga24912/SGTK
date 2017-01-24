//
// Created by olga on 10.10.16.
//

#ifndef SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H
#define SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

#include <bits/stdc++.h>

using namespace std;

/*
 * system call of other programs
 */
class SystemAlignmentTools {
public:
    /*
     * alignment RNA reads from rnaFile to refFile
     * result will be in resFile.
     */
    void alignmentRNA(string refFileName, string rnaFileName, string resFileName);
    static void alignmentREF(string refFileName, string queryFileName);
};


#endif //SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

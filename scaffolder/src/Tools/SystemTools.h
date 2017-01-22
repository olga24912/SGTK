//
// Created by olga on 10.10.16.
//

#ifndef SCAFFOLDER_SYSTEMTOOLS_H
#define SCAFFOLDER_SYSTEMTOOLS_H

#include <bits/stdc++.h>

using namespace std;

/*
 * system call of other programs
 */
class SystemTools {
public:
    /*
     * alignment RNA reads from rnaFile to refFile
     * result will be in resFile.
     */
    void alignmentRNA(string refFileName, string rnaFileName, string resFileName);
    static void alignmentREF(string refFileName, string queryFileName);
    static void showDotFile(string fileName);
};


#endif //SCAFFOLDER_SYSTEMTOOLS_H

//
// Created by olga on 10.10.16.
//

#ifndef SCAFFOLDER_WORKWITHOTHERTOOLS_H
#define SCAFFOLDER_WORKWITHOTHERTOOLS_H

#include <bits/stdc++.h>

using namespace std;

/*
 * system call of other programs
 */
class WorkWithOtherTools {
public:
    /*
     * alignment RNA reads from rnaFile to refFile
     * result will be in resFile.
     */
    void alignmentRNA(string refFileName, string rnaFileName, string resFileName);
};


#endif //SCAFFOLDER_WORKWITHOTHERTOOLS_H

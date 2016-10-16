//
// Created by olga on 16.10.16.
//

#ifndef SCAFFOLDER_SPLITREADSTEST_H
#define SCAFFOLDER_SPLITREADSTEST_H

#include <bits/stdc++.h>
using namespace std;

class SplitReadsTest {
private:
    string refOUTFileName;
    string refINFileName;
    string readsFileName;

    int readLen;

    string smallRef(string ref, int exonLen, int intronLen, int spos);
public:
    void genTest(string fi, string fo, string rf, int readLen);
};


#endif //SCAFFOLDER_SPLITREADSTEST_H

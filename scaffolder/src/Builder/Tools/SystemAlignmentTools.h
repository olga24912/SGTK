#ifndef SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H
#define SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

#include <bits/stdc++.h>

//system call of other programs
class SystemAlignmentTools {
public:
    // alignment RNA reads from rnaFile to refFile
    // result will be in resFile.
    void alignmentRNA(std::string refFileName, std::string rnaFileName, std::string resFileName, std::string path=".");
    static void alignmentREF(std::string refFileName, std::string queryFileName);
};


#endif //SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

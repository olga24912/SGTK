#ifndef SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H
#define SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

#include <string>

//system call of other programs for alignment
class SystemAlignmentTools {
public:
    // alignment RNA reads from rnaFile to refFile
    // result will be in resFile.
    static void alignmentRNA(std::string refFileName, std::string rnaFileName, std::string resFileName, std::string path=".");
    // alignment queryFile on RefFile result will be in "out.coords" file.
    static void alignmentREF(std::string refFileName, std::string queryFileName);
};


#endif //SCAFFOLDER_SYSTEMALIGNMENTTOOLS_H

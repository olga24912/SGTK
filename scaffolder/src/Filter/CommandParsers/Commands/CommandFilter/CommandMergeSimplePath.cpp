#include "CommandMergeSimplePath.h"

void CommandMergeSimplePath::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string contigFile;
    std::string outFile;

    ss >> contigFile >> outFile;

    ScafSimplePath ssp;
    ssp.evaluate(filter, contigFile, outFile);

}

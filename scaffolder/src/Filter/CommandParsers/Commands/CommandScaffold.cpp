#include <Filter/Scaffolder/ScaffolderPipeline.h>
#include "CommandScaffold.h"

void CommandScaffold::execute(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string contigFile;
    std::string outFile;

    ss >> contigFile >> outFile;

    ScaffolderPipeline sp;
    sp.evaluate(filter, contigFile, outFile);

}

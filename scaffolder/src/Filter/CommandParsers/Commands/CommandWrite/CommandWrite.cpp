#include "CommandWrite.h"

void CommandWrite::execute(std::string argv, State& state, Filter *filter) {
    std::string fileName;
    std::stringstream ss(argv);
    ss >> fileName;

    writeGraph(argv, state, filter);

    SystemTools::showDotFile(fileName);
}

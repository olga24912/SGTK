#include "CommandWriteFull.h"

void CommandWriteFull::writeGraph(std::string argv, State& state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    ss >> fileName;

    INFO("write full graph fileName=" << fileName);

    WriteFullGraph writer(fileName, filter, state.validator, state.dotWriterBuilder);
    writer.write();
}

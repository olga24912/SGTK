#include "CommandWriteBigComp.h"

void CommandWriteBigComp::writeGraph(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int size;
    ss >> fileName >> size;
    WriteBigComponent writer(fileName, size, filter, state.validator);
    writer.write();
}

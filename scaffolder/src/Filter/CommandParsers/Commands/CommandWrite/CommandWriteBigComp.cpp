#include "CommandWriteBigComp.h"

void CommandWriteBigComp::writeGraph(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int size;
    ss >> fileName >> size;

    INFO("Write Big Comp fileName=" << fileName << " size=" << size);

    WriteBigComponent writer(fileName, size, filter, state.validator, state.dotWriterBuilder);
    writer.write();
}

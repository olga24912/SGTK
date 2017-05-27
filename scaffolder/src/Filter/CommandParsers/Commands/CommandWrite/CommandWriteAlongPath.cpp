#include "CommandWriteAlongPath.h"

// writeAlongPath <prefixFileName> <libNum> <dist> <minRefPathSize>
void CommandWriteAlongPath::writeGraph(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    ss >> fileName;
    int libNum;
    int dist;
    int minSize;
    ss >> libNum >> dist >> minSize;

    WriteAlongPath writer(fileName, libNum, dist, minSize, filter, state.validator, state.dotWriterBuilder);

    writer.write();
}

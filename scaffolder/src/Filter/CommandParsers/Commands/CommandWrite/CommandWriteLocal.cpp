#include "CommandWriteLocal.h"

void CommandWriteLocal::writeGraph(std::string argv, State& state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int v;
    int dist;
    ss >> fileName >> v >> dist;

    INFO("write local graph fileName=" << fileName << " v=" << v << " dist=" << dist);

    WriteLocal writer(v, dist, fileName, filter, state.validator, state.dotWriterBuilder);
    writer.write();

    state.name = State::LOCAL;
    state.fileName = fileName;
    state.dist = dist;
}

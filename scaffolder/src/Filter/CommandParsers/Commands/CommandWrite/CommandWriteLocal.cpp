#include "CommandWriteLocal.h"

void CommandWriteLocal::writeGraph(std::string argv, State& state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int v;
    int dist;
    ss >> fileName >> v >> dist;

    std::cerr << "command write local" << std::endl;

    WriteLocal writer(v, dist, fileName, filter, state.validator, state.maxVert, state.maxEdge);
    writer.write();

    state.name = State::LOCAL;
    state.fileName = fileName;
    state.dist = dist;
}

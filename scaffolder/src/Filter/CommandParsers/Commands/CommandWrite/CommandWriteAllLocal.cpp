#include "CommandWriteAllLocal.h"

void CommandWriteAllLocal::writeGraph(std::string argv, State &state, Filter *filter) {
    std::stringstream ss(argv);
    std::string fileName;
    int dist;
    ss >> fileName >> dist;

    INFO("write all local fileName=" << fileName << " dist=" << dist);

    for (int v = 0; v < filter->getVertexCount(); ++v) {
        string name = "";
        int x = v;
        while (x > 0) {
            name += '0' + (x%10);
            x /= 10;
        }
        reverse(name.begin(), name.end());
        name = fileName + name;

        WriteLocal writer(v, dist, name, filter, state.validator, state.dotWriterBuilder);
        writer.write();
    }

    state.name = State::LOCAL;
    state.fileName = fileName;
    state.dist = dist;
}
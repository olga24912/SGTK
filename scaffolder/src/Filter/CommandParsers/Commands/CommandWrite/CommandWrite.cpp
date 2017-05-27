#include "CommandWrite.h"

void CommandWrite::execute(std::string argv, State& state, Filter *filter) {
    std::string fileName;
    std::stringstream ss(argv);
    ss >> fileName;

    state.dotWriterBuilder.setFilter(filter);
    state.dotWriterBuilder.setValidator(state.validator);
    state.dotWriterBuilder.setMaxVert(state.maxVert);
    state.dotWriterBuilder.setMaxEdge(state.maxEdge);

    writeGraph(argv, state, filter);

    SystemTools::showDotFile(fileName);
}

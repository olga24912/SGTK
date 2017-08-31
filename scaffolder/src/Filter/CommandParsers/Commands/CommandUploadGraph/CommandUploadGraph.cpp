#include <Filter/CommandParsers/State.h>
#include "CommandUploadGraph.h"

namespace filter {
    namespace commands {
        void CommandUploadGraph::execute(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            ss >> fileName;

            INFO("uploadGraph fileName=" << fileName);
            graph = ContigGraph::read(fileName);
        }
    }
}

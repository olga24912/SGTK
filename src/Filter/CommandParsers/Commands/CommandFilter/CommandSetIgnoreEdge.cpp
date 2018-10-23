#include <Filter/ContigGraph/ContigGraph.h>
#include "CommandSetIgnoreEdge.h"

namespace filter {
    namespace commands {
        void CommandSetIgnoreEdge::execute(std::string argv, State &state, ContigGraph &graph) {
            INFO("set graph set ignore edge");
            std::stringstream ss;
            int e;
            ss >> e;
            graph.delEdge(e);
        }
    }
}

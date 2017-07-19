#include <Filter/CommandParsers/State.h>
#include "CommandSetIgnore.h"

namespace filter {
    namespace commands {
        void CommandSetIgnore::execute(std::string argv, State &state, ContigGraph &graph) {
            INFO("set graph set ignore");
            int vb, ve;
            std::stringstream ss(argv);
            ss >> vb >> ve;
            for (int v = vb; v < ve; ++v) {
                graph.delVertex(v);
            }
        }
    }
}

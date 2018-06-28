#include <Filter/CommandParsers/State.h>
#include "CommandMinContig.h"

namespace filter {
    namespace commands {
        void CommandMinContig::execute(std::string argv, State &state, ContigGraph &graph) {
            INFO("set ContigGraph min contig");
            int minLen;
            std::stringstream ss(argv);
            ss >> minLen;

            std::vector<int> vert = graph.getVertexList();
            for (int v : vert) {
                if (graph.getTargetLen(v) < minLen) {
                    graph.delVertex(v);
                }
            }
        }
    }
}
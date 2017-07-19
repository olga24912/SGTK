#include "CommandMergeLib.h"

namespace filter {
    namespace commands {
        void CommandMergeLib::execute(std::string argv, State &state, ContigGraph &graph) {
            INFO("set ContigGraph merge lib");
            int l1, l2;
            std::stringstream ss(argv);
            ss >> l1 >> l2;
            std::string name;
            ss >> name;
            graph.mergeLib(l1, l2, name, 1, 1);
        }
    }
}
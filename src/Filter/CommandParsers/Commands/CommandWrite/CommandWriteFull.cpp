#include "CommandWriteFull.h"

namespace filter {
    namespace commands {
        void CommandWriteFull::writeGraph(std::string argv, State &state, ContigGraph &graph) {
                std::stringstream ss(argv);
                std::string fileName;
                ss >> fileName;

                INFO("write full graph fileName=" << fileName);

                writers::WriteFullGraph writer(fileName, &graph, state.validator, state.dotWriterBuilder);
                writer.write();
        }
    }
}
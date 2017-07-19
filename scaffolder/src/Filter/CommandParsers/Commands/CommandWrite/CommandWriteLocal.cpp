#include "CommandWriteLocal.h"

namespace filter {
    namespace commands {
        void CommandWriteLocal::writeGraph(std::string argv, State &state, ContigGraph &graph) {
                std::stringstream ss(argv);
                std::string fileName;
                int v;
                int dist;
                ss >> fileName >> v >> dist;

                INFO("write local graph fileName=" << fileName << " v=" << v << " dist=" << dist);

                writers::WriteLocal writer(v, dist, fileName, &graph, state.validator, state.dotWriterBuilder);
                writer.write();

                state.name = State::LOCAL;
                state.fileName = fileName;
                state.dist = dist;
        }
    }
}
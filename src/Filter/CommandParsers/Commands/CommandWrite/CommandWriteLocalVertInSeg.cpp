#include "CommandWriteLocalVertInSeg.h"

namespace filter {
    namespace commands {
        void CommandWriteLocalVertInSeg::writeGraph(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            int vb, ve;
            int dist;
            ss >> fileName >> vb >> ve >> dist;

            INFO("write local vert in seg fileName=" << fileName << " vb=" << vb << " ve=" << ve << " dist=" << dist);

            for (int v = vb; v <= ve; ++v) {
                std::string name = "";
                int x = v;
                while (x > 0) {
                    name += '0' + (x % 10);
                    x /= 10;
                }
                std::reverse(name.begin(), name.end());
                name = fileName + name;

                writers::WriteLocal writer(v, dist, name, &graph, state.validator, state.dotWriterBuilder);
                writer.write();
            }

            state.name = State::LOCAL;
            state.fileName = fileName;
            state.dist = dist;
        }
    }
}
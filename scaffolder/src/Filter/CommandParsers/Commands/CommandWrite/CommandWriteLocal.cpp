#include "CommandWriteLocal.h"

namespace filter {
    namespace commands {
        void CommandWriteLocal::writeGraph(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            std::vector<int> v;
            int dist;
            ss >> fileName >> dist;
            int w;
            while (ss >> w) {
                v.push_back(w);
            }

            INFO("write local graph fileName=" << fileName << " len(v)=" << v.size() << " dist=" << dist);

            for (int w : v) {
                std::stringstream sw;
                sw << fileName.c_str();
                sw << "_" << w << "_";
                writers::WriteLocal writer(w, dist, sw.str(), &graph, state.validator, state.dotWriterBuilder);
                writer.write();
            }

            state.name = State::LOCAL;
            state.fileName = fileName;
            state.dist = dist;
        }
    }
}
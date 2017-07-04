#include <Filter/CommandParsers/State.h>
#include "CommandUploadGraph.h"

namespace filter {
    namespace commands {
        void CommandUploadGraph::execute(std::string argv, State &state, Filter *filter) {
            std::stringstream ss(argv);
            std::string fileName;
            ss >> fileName;

            INFO("uploadGraph fileName=" << fileName);
            filter->processQuery(Query(Query::UPLOAD_GRAPH, fileName));
        }
    }
}

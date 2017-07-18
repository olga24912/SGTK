#include <Filter/CommandParsers/State.h>
#include "CommandMinEdgeWeight.h"

namespace filter {
    namespace commands {
        void CommandMinEdgeWeight::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set ContigGraph min Edge weight");
            filter->processQuery(Query(Query::MIN_EDGE_WEIGHT, argv));
        }
    }
}

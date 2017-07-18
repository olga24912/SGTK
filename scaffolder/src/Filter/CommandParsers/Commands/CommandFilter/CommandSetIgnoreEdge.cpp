#include "CommandSetIgnoreEdge.h"

namespace filter {
    namespace commands {

        void CommandSetIgnoreEdge::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set graph set ignore edge");
            filter->processQuery(Query(Query::SET_IGNORE_EDGE, argv));
        }
    }
}

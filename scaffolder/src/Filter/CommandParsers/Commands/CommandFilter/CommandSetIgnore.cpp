#include <Filter/CommandParsers/State.h>
#include "CommandSetIgnore.h"

namespace filter {
    namespace commands {
        void CommandSetIgnore::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set graph set ignore");
            filter->processQuery(Query(Query::SET_IGNORE, argv));
        }
    }
}

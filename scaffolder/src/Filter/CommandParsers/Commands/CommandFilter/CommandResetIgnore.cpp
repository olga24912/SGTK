#include <Filter/CommandParsers/State.h>
#include "CommandResetIgnore.h"

namespace filter {
    namespace commands {
        void CommandResetIgnore::execute(std::string argv, State &state, Filter *filter) {
            INFO("set Filter reset ignore");
            filter->processQuery(Query(Query::RESET_IGNORE, argv));
        }
    }
}

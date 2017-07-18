#include <Filter/CommandParsers/State.h>
#include "CommandResetIgnore.h"

namespace filter {
    namespace commands {
        void CommandResetIgnore::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set ContigGraph reset ignore");
            filter->processQuery(Query(Query::RESET_IGNORE, argv));
        }
    }
}

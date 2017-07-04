#include <Filter/CommandParsers/State.h>
#include "CommandPrint.h"

namespace filter {
    namespace commands {
        void CommandPrint::execute(std::string argv, State &state, Filter *filter) {
            INFO("print graph");
            filter->write(argv);
        }
    }
}
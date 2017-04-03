#include <Filter/CommandParsers/State.h>
#include "CommandPrint.h"

void CommandPrint::execute(std::string argv, State &state, Filter *filter) {
    filter->write(argv);
}

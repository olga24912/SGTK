//
// Created by olga on 26.01.17.
//

#include <Filter/CommandParsers/State.h>
#include "CommandResetIgnore.h"

void CommandResetIgnore::execute(std::string argv, State& state, Filter *filter) {
    filter->processQuery(Query(Query::RESET_IGNORE, argv));
}

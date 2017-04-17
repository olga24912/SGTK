//
// Created by olga on 26.01.17.
//

#include <Filter/CommandParsers/State.h>
#include "CommandSetIgnore.h"

void CommandSetIgnore::execute(std::string argv, State& state, Filter *filter) {
    filter->processQuery(Query(Query::SET_IGNORE, argv));
}

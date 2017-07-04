#include <Filter/CommandParsers/State.h>
#include "CommandSetIgnore.h"

void CommandSetIgnore::execute(std::string argv, State& state, Filter *filter) {
    INFO("set filter set ignore");
    filter->processQuery(Query(Query::SET_IGNORE, argv));
}

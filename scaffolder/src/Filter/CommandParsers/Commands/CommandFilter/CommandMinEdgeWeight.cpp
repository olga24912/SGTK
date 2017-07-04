#include <Filter/CommandParsers/State.h>
#include "CommandMinEdgeWeight.h"

void CommandMinEdgeWeight::execute(std::string argv, State& state, Filter *filter) {
    INFO("set Filter min Edge weight");
    filter->processQuery(Query(Query::MIN_EDGE_WEIGHT, argv));
}

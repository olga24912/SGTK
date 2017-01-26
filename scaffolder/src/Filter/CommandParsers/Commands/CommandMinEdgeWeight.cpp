//
// Created by olga on 26.01.17.
//

#include <Filter/CommandParsers/State.h>
#include "CommandMinEdgeWeight.h"

void CommandMinEdgeWeight::execute(std::string argv, State& state, Filter *filter) {
    filter->processQuery(Query(Query::MIN_EDGE_WEIGHT, argv));
}

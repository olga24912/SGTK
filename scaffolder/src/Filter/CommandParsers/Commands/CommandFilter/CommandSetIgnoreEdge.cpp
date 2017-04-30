#include "CommandSetIgnoreEdge.h"

void CommandSetIgnoreEdge::execute(std::string argv, State &state, Filter *filter) {
    filter->processQuery(Query(Query::SET_IGNORE_EDGE, argv));
}

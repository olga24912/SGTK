#include "CommandSetIgnoreEdge.h"

void CommandSetIgnoreEdge::execute(std::string argv, State &state, Filter *filter) {
    INFO("set filter set ignore edge");
    filter->processQuery(Query(Query::SET_IGNORE_EDGE, argv));
}

#include "CommandMergeLib.h"


void CommandMergeLib::execute(std::string argv, State& state, Filter *filter) {
    INFO("set Filter merge lib");
    filter->processQuery(Query(Query::MERGE_LIB, argv));
}
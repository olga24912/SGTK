#include "CommandMergeLib.h"


void CommandMergeLib::execute(std::string argv, State& state, Filter *filter) {
    filter->processQuery(Query(Query::MERGE_LIB, argv));
}
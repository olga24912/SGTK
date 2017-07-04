#include <Filter/CommandParsers/State.h>
#include "CommandMinContig.h"

void CommandMinContig::execute(std::string argv, State& state, Filter *filter) {
    INFO("set Filter min contig");
    filter->processQuery(Query(Query::MIN_CONTIG_LEN, argv));
}

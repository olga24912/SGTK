//
// Created by olga on 26.01.17.
//

#include <Filter/CommandParsers/State.h>
#include "CommandMinContig.h"

void CommandMinContig::execute(std::string argv, State state, Filter *filter) {
    filter->processQuery(Query(Query::MIN_CONTIG_LEN, argv));
}

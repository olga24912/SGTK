#include <Filter/CommandParsers/State.h>
#include "CommandMinContig.h"

namespace filter {
    namespace commands {
        void CommandMinContig::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set ContigGraph min contig");
            filter->processQuery(Query(Query::MIN_CONTIG_LEN, argv));
        }
    }
}
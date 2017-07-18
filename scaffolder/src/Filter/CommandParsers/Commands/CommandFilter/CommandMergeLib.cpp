#include "CommandMergeLib.h"

namespace filter {
    namespace commands {
        void CommandMergeLib::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("set ContigGraph merge lib");
            filter->processQuery(Query(Query::MERGE_LIB, argv));
        }
    }
}
#include "CommandMergeLib.h"

namespace filter {
    namespace commands {
        void CommandMergeLib::execute(std::string argv, State &state, Filter *filter) {
            INFO("set Filter merge lib");
            filter->processQuery(Query(Query::MERGE_LIB, argv));
        }
    }
}
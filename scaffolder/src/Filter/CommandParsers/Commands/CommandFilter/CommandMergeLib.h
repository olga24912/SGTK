#ifndef SCAFFOLDER_COMMANDMERGELIB_H
#define SCAFFOLDER_COMMANDMERGELIB_H


#include <Filter/Filters/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        class CommandMergeLib : public Command {
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDMERGELIB_H

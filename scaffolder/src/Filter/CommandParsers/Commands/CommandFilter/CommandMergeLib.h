#ifndef SCAFFOLDER_COMMANDMERGELIB_H
#define SCAFFOLDER_COMMANDMERGELIB_H


#include <Filter/Filters/Filter.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        class CommandMergeLib : public Command {
            void execute(std::string argv, State &state, Filter *filter) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDMERGELIB_H

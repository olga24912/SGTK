#ifndef SCAFFOLDER_COMMANDMERGELIB_H
#define SCAFFOLDER_COMMANDMERGELIB_H


#include <Filter/ContigGraph/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        /*
         * mergeLib <libNum1> <libNum2> <newLibName>
         *
         * Merge two source with id <libNum1> and id <libNum2>. New lib will have id <libNum1> and name <newLibName>
         */

        class CommandMergeLib : public Command {
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDMERGELIB_H

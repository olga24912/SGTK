#ifndef SCAFFOLDER_COMMANDSETIGNORE_H
#define SCAFFOLDER_COMMANDSETIGNORE_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        /*
         * setIgnore <vertexIdStart> <vertexIdFinish>
         *
         * Delete all node with ID between <vertexIdStart> and <vertexIdFinish>
         */

        class CommandSetIgnore : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}


#endif //SCAFFOLDER_COMMANDSETIGNORE_H

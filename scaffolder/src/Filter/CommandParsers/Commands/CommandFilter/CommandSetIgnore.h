#ifndef SCAFFOLDER_COMMANDSETIGNORE_H
#define SCAFFOLDER_COMMANDSETIGNORE_H

#include <Filter/Filters/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {

        class CommandSetIgnore : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}


#endif //SCAFFOLDER_COMMANDSETIGNORE_H

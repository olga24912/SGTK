#ifndef SCAFFOLDER_COMMANDRESETIGNORE_H
#define SCAFFOLDER_COMMANDRESETIGNORE_H

#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        class CommandResetIgnore : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}


#endif //SCAFFOLDER_COMMANDRESETIGNORE_H

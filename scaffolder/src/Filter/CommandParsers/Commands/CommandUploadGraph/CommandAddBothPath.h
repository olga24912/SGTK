#ifndef SCAFFOLDER_COMMANDADDBOTHPATH_H
#define SCAFFOLDER_COMMANDADDBOTHPATH_H


#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandAddBothPath : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;

        };
    }
}


#endif //SCAFFOLDER_COMMANDADDBOTHPATH_H

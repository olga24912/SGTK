#ifndef SCAFFOLDER_COMMANDSETIGNOREEDGE_H
#define SCAFFOLDER_COMMANDSETIGNOREEDGE_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandSetIgnoreEdge : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}


#endif //SCAFFOLDER_COMMANDSETIGNOREEDGE_H

#ifndef SCAFFOLDER_COMMANDWRONGRIGHTSTATISTIC_H
#define SCAFFOLDER_COMMANDWRONGRIGHTSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandWrongRightStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;

        };
    }
}


#endif //SCAFFOLDER_COMMANDWRONGRIGHTSTATISTIC_H

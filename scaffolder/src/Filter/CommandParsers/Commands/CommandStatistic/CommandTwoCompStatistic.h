#ifndef SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H
#define SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {

        class CommandTwoCompStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDTWOCOMPSTATISTIC_H

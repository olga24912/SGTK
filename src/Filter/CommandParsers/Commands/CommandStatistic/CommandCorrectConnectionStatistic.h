#ifndef SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H
#define SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H


#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {

        class CommandCorrectConnectionStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDCORRECTCONNECTIONSTATISTIC_H

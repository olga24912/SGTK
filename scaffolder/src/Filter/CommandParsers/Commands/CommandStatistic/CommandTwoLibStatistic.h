#ifndef SCAFFOLDER_COMMANDTWOLIBSTATISTIC_H
#define SCAFFOLDER_COMMANDTWOLIBSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandTwoLibStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}
#endif //SCAFFOLDER_COMMANDTWOLIBSTATISTIC_H

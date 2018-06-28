#ifndef SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H
#define SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H

#include <Filter/Statistics/Statistic.h>
#include <Filter/CommandParsers/State.h>
#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandWeightStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDWEIGHTSTATISTIC_H

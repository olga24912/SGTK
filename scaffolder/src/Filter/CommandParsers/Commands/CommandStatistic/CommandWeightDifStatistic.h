#ifndef SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H
#define SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {

        class CommandWeightDifStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}
#endif //SCAFFOLDER_COMMANDWEIGHTDIFSTATISTIC_H

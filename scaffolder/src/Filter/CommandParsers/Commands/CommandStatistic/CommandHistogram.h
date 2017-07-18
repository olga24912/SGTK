#ifndef SCAFFOLDER_COMMANDHISTOGRAM_H
#define SCAFFOLDER_COMMANDHISTOGRAM_H

#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {

        class CommandHistogram : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}
#endif //SCAFFOLDER_COMMANDHISTOGRAM_H

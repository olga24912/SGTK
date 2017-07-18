#ifndef SCAFFOLDER_COMMANDPRINT_H
#define SCAFFOLDER_COMMANDPRINT_H

#include <Filter/Filters/ContigGraph.h>
#include "Command.h"

namespace filter {
    namespace commands {
        class CommandPrint : public Command {

        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDPRINT_H

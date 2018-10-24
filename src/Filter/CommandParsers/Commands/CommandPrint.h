#ifndef SCAFFOLDER_COMMANDPRINT_H
#define SCAFFOLDER_COMMANDPRINT_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "Command.h"

namespace filter {
    namespace commands {
       /*
        * print <fileName>
        *
        * Serialize graph to <fileName> in .gr format
        */

        class CommandPrint : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDPRINT_H

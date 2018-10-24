#ifndef SCAFFOLDER_COMMANDMERGESIMPLEPATH_H
#define SCAFFOLDER_COMMANDMERGESIMPLEPATH_H

#include <Filter/CommandParsers/State.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        /*
         * mergeSimplePath <contigsFileName> <outFileName>
         *
         * Build scaffolds. Scaffolds will be output to <outFileName> and build from contigs from <contigsFileName>
         */

        class CommandScaffold : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}
#endif //SCAFFOLDER_COMMANDMERGESIMPLEPATH_H

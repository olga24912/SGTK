#ifndef SCAFFOLDER_COMMANDMINCONTIG_H
#define SCAFFOLDER_COMMANDMINCONTIG_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        /*
         * minContig <len>
         *
         * Delete all contig with length less then <len>
         */

        class CommandMinContig : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDMINCONTIG_H

#ifndef SCAFFOLDER_COMMANDMINEDGEWEIGHT_H
#define SCAFFOLDER_COMMANDMINEDGEWEIGHT_H

#include <Filter/ContigGraph/ContigGraph.h>
#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        /*
         * minEdgeW <libNum> <weight>
         *
         * Delete all edges from source with ID <libNum>, with wieght less then <weight>
         */

        class CommandMinEdgeWeight : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDMINEDGEWEIGHT_H

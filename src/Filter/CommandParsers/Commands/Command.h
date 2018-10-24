#ifndef SCAFFOLDER_COMMAND_H
#define SCAFFOLDER_COMMAND_H

#include <string>
#include "Filter/CommandParsers/State.h"
#include <sstream>
#include <Filter/ContigGraph/ContigGraph.h>

namespace filter {
    namespace commands {
        using namespace contig_graph;
        /*
         * Interface for commands which will be apply
         */
        class Command {
        public:
            // apply command
            virtual void execute(std::string argv, State &state, ContigGraph &graph) = 0;

        protected:
            DECL_LOGGER("Command");
        };
    }
}

#endif //SCAFFOLDER_COMMAND_H

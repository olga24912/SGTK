#ifndef SCAFFOLDER_COMMANDUPLOADGRAPH_H
#define SCAFFOLDER_COMMANDUPLOADGRAPH_H

#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
       /*
        * uploadGraph <filename>
        *
        * Upload graph from <filename> in .gr format
        */

       class CommandUploadGraph : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDUPLOADGRAPH_H

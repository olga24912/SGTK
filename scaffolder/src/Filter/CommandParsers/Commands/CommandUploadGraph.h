#ifndef SCAFFOLDER_COMMANDUPLOADGRAPH_H
#define SCAFFOLDER_COMMANDUPLOADGRAPH_H

#include <Filter/Filters/ContigGraph.h>
#include "Command.h"

namespace filter {
    namespace commands {
        //uploadGraph <filename>
        class CommandUploadGraph : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDUPLOADGRAPH_H

#ifndef SCAFFOLDER_COMMANDWRITEALONGPATH_H
#define SCAFFOLDER_COMMANDWRITEALONGPATH_H

#include "CommandWrite.h"

namespace filter {
    namespace commands {
        class CommandWriteAlongPath : public CommandWrite {
        public:
            void writeGraph(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDWRITEALONGPATH_H

#ifndef SCAFFOLDER_COMMANDWRITELOCAL_H
#define SCAFFOLDER_COMMANDWRITELOCAL_H

#include "CommandWrite.h"
#include <Filter/Writers/WriteLocal.h>

namespace filter {
    namespace commands {
        class CommandWriteLocal : public CommandWrite {
        public:
            void writeGraph(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDWRITELOCAL_H

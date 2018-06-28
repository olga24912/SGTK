#ifndef SCAFFOLDER_COMMANDWRITEBIGCOMP_H
#define SCAFFOLDER_COMMANDWRITEBIGCOMP_H

#include "CommandWrite.h"
#include <Filter/Writers/WriteBigComponent.h>

namespace filter {
    namespace commands {
        class CommandWriteBigComp : public CommandWrite {
        public:
            void writeGraph(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}

#endif //SCAFFOLDER_COMMANDWRITEBIGCOMP_H

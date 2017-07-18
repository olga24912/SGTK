#ifndef SCAFFOLDER_COMMANDWRITE_H
#define SCAFFOLDER_COMMANDWRITE_H

#include <Filter/Filters/ContigGraph.h>
#include <Filter/CommandParsers/State.h>
#include <Filter/Tools/SystemTools.h>
#include <Filter/CommandParsers/Commands/Command.h>

namespace filter {
    namespace commands {
        class CommandWrite : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph *filter) override;

            virtual void writeGraph(std::string argv, State &state, ContigGraph *filter) = 0;
        };
    }
}

#endif //SCAFFOLDER_COMMANDWRITE_H

#ifndef SCAFFOLDER_COMMANDADDBOTHPATH_H
#define SCAFFOLDER_COMMANDADDBOTHPATH_H


#include <Filter/CommandParsers/Commands/Command.h>
#include <Filter/ContigGraph/ContigGraph.h>

namespace filter {
    namespace commands {
        /*
         * addBothPath <bothPathFile> <newLibName> <colorLib>
         *
         * Add scaffolds to graph from <bothPathFile>
         */

        class CommandAddBothPath : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;

        };
    }
}


#endif //SCAFFOLDER_COMMANDADDBOTHPATH_H

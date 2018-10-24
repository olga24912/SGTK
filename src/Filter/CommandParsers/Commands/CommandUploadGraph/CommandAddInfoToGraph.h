#ifndef SCAFFOLDER_COMMANDADDINFOTOGRAPH_H
#define SCAFFOLDER_COMMANDADDINFOTOGRAPH_H


#include <Filter/CommandParsers/Commands/Command.h>
#include <Filter/ContigGraph/ContigGraph.h>

namespace filter {
    namespace commands {
        /*
         * addInfoToGraph <infoFile> <newLibName> <colorLib>
         *
         * Add scaffolds to graph from <infoFile>
         */

        class CommandAddInfoToGraph : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;
        };
    }
}


#endif //SCAFFOLDER_COMMANDADDINFOTOGRAPH_H

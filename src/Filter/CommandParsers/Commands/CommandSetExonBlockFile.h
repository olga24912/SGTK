#ifndef SCAFFOLDER_COMMANDSETEXONBLOCKFILE_H
#define SCAFFOLDER_COMMANDSETEXONBLOCKFILE_H

#include "Command.h"

namespace filter {
    namespace commands {
       /*
        * setExonBlock <crdFile>
        *
        * Set file with annotation in crd format
        */

        class CommandSetExonBlockFile : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                std::string fileName;
                std::stringstream ss(argv);
                ss >> fileName;
                
                INFO("set exon block file");
                graph.addExonBlock(fileName);
            }
        };
    }
}


#endif //SCAFFOLDER_COMMANDSETEXONBLOCKFILE_H

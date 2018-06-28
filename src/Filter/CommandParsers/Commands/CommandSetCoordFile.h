#ifndef SCAFFOLDER_COMMANDSETCOORDFILE_H
#define SCAFFOLDER_COMMANDSETCOORDFILE_H

#include "Command.h"

namespace filter {
    namespace commands {
        class CommandSetCoordFile : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                std::string coordFile;
                std::stringstream ss(argv);
                ss >> coordFile;
                INFO("set coord File");
                state.coordFile = coordFile;
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSETCOORDFILE_H

#ifndef SCAFFOLDER_COMMANDUPLOADSCAFFOLDSGRAPH_H
#define SCAFFOLDER_COMMANDUPLOADSCAFFOLDSGRAPH_H

#include "Filter/CommandParsers/Commands/Command.h"

namespace filter {
    namespace commands {
        class CommandUploadScaffoldsGraph : public Command {
        private:
            std::string getContigName(std::string s);
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override;

            void addVertexis(ContigGraph &graph, int libNum, std::string contigFileName);
        };
    }
}


#endif //SCAFFOLDER_COMMANDUPLOADSCAFFOLDSGRAPH_H

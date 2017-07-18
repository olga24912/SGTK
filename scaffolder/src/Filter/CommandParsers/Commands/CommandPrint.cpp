#include <Filter/CommandParsers/State.h>
#include "CommandPrint.h"

namespace filter {
    namespace commands {
        void CommandPrint::execute(std::string argv, State &state, ContigGraph *filter) {
            INFO("print graph");
            std::string fileName;
            std::stringstream ss(argv);
            ss >> fileName;
            filter->write(fileName);
        }
    }
}
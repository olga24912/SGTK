#include <Filter/Scaffolder/ScaffolderPipeline.h>
#include "CommandScaffold.h"

namespace filter {
    namespace commands {
        void CommandScaffold::execute(std::string argv, State &state, ContigGraph *filter) {
            std::stringstream ss(argv);
            std::string contigFile;
            std::string outFile;

            ss >> contigFile >> outFile;

            INFO("create scaffolds contigFile=" << contigFile << " outFile=" << outFile);

            scaffolder::ScaffolderPipeline sp;
            sp.evaluate(filter, contigFile, outFile);
        }
    }
}

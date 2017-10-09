#include <Filter/Scaffolder/ScaffolderPipeline.h>
#include "CommandScaffold.h"

namespace filter {
    namespace commands {
        void CommandScaffold::execute(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string contigFile;
            std::string outFile;

            ss >> contigFile >> outFile;

            INFO("create scaffolds contigFile=" << contigFile << " outFile=" << outFile);

            scaffolder::ScaffolderPipeline sp;

            alig_info::InfoAboutContigsAlig alig;
            alig.parseCoordFile(&graph, state.coordFile);

            sp.setAlig(&alig);
            sp.evaluate(&graph, contigFile, outFile, state.bamFiles);
        }
    }
}

#ifndef SCAFFOLDER_COMMANDSETBAMFILE_H
#define SCAFFOLDER_COMMANDSETBAMFILE_H

#include "Command.h"

namespace filter {
    namespace commands {
        class CommandSetBamFile : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                std::string bam1, bam2, bai1, bai2;
                int lib;
                std::stringstream ss(argv);
                State::BamFiles bamFiles;
                ss >> bamFiles.bam1 >> bamFiles.bam2 >> bamFiles.bai1 >> bamFiles.bai2 >> bamFiles.lib;
                INFO("set bam files");
                state.bamFiles.push_back(bamFiles);
            }
        };
    }
}


#endif //SCAFFOLDER_COMMANDSETBAMFILE_H

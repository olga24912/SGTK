#ifndef SCAFFOLDER_STATE_H
#define SCAFFOLDER_STATE_H

#include <string>
#include <vector>

namespace filter {
    namespace commands {
        struct State {
            struct BamFiles {
                std::string bam1;
                std::string bam2;
                std::string bai1;
                std::string bai2;
                int lib;
            };

            std::string fileName = "";
            std::string coordFile = "";

            std::vector<BamFiles> bamFiles;
        };
    }
}
#endif //SCAFFOLDER_STATE_H

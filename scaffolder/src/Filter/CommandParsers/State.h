#ifndef SCAFFOLDER_STATE_H
#define SCAFFOLDER_STATE_H

#include <string>
#include <Filter/Writers/FileValidator/FileValidator.h>
#include <Filter/Writers/DotWriter/DotWriter.h>
#include <Filter/Writers/DotWriter/DotWriterBuilder.h>

namespace filter {
    namespace commands {
        struct State {
            enum StateName {
                DEF, LOCAL
            };

            struct BamFiles {
                std::string bam1;
                std::string bam2;
                std::string bai1;
                std::string bai2;
                int lib;
            };

            StateName name = StateName::DEF;
            std::string fileName = "";
            int dist = 0;
            writers::FileValidator *validator = new writers::FileValidator;
            writers::DotWriterBuilder *dotWriterBuilder = new writers::DotWriterBuilder();
            int maxEdge = 40;
            int maxVert = 20;
            std::string coordFile = "";

            std::vector<BamFiles> bamFiles;

            ~State() {
                delete validator;
            }
        };
    }
}
#endif //SCAFFOLDER_STATE_H

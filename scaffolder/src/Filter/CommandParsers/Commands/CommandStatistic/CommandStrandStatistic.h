#ifndef SCAFFOLDER_COMMANDSTRANDSTATISTIC_H
#define SCAFFOLDER_COMMANDSTRANDSTATISTIC_H


#include <Filter/CommandParsers/Commands/Command.h>
#include <Filter/Statistics/StrandStatistic.h>

namespace filter {
    namespace commands {
        class CommandStrandStatistic : public Command {
        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                std::stringstream ss(argv);
                std::string gffFile;
                std::string coordFile;
                ss >> gffFile >> coordFile;

                INFO("STRAND STATISTIC");

                statistics::StrandStatistic cs;
                cs.calculateStatistic(&graph, coordFile, gffFile);
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSTRANDSTATISTIC_H

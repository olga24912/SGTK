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
                std::string gffFileR, gffFileC, out;
                std::string coordFile;
                ss >> gffFileR >> coordFile >> gffFileC >> out;

                INFO("STRAND STATISTIC");

                alig_info::InfoAboutContigsAlig aligInfo;
                aligInfo.parseCoordFile(&graph, coordFile);
                statistics::StrandStatistic cs(&graph, aligInfo, gffFileR, gffFileC, out);
                cs.calculateStatistic();
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSTRANDSTATISTIC_H

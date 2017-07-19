#include <Filter/Statistics/CntRightConnectionStatistic.h>
#include "CommandCorrectConnectionStatistic.h"

namespace filter {
    namespace commands {
        void CommandCorrectConnectionStatistic::execute(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            int libNum;
            ss >> fileName >> libNum;

            INFO("Correct connection statistic fileName=" << fileName << " libNum=" << libNum);

            statistics::CntRightConnectionStatistic cs;
            cs.calculateStatistic(&graph, fileName, libNum);
        }
    }
}
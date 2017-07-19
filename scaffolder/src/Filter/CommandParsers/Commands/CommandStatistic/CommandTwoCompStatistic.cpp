#include <Filter/Statistics/TwoCompStatistic.h>
#include "CommandTwoCompStatistic.h"


namespace filter {
    namespace commands {
        void CommandTwoCompStatistic::execute(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            int libNum1;
            int val1;
            int libNum2;
            int val2;
            ss >> fileName >> libNum1 >> val1 >> libNum2 >> val2;

            INFO("TwoComp statistic libNum1=" << libNum1 << " val1=" << val1 << " libNum2=" << libNum2 << " val2="
                                              << val2);

            statistics::TwoCompStatistic cs;
            cs.calculateStatistic(&graph, fileName, libNum1, val1, libNum2, val2);
        }
    }
}
#include <Filter/Statistics/DifStatistic.h>
#include "CommandDifStatistic.h"

namespace filter {
    namespace commands {

        void CommandDifStatistic::execute(std::string argv, State &state, ContigGraph &graph) {
                std::stringstream ss(argv);
                std::string fileName;
                int libNum;
                int step;
                int mxVal;
                ss >> fileName >> libNum >> step >> mxVal;

                INFO("Dif statistic fileName=" << fileName << " libNum=" << libNum << " step=" << step << " mxVal="
                                               << mxVal);

                statistics::DifStatistic cs;
                cs.calculateStatistic(&graph, fileName, libNum, step, mxVal);
        }
    }
}
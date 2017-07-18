#include <Filter/Statistics/WeightStatistic.h>
#include "CommandWeightStatistic.h"

namespace filter {
    namespace commands {
        void CommandWeightStatistic::execute(std::string argv, State &state, ContigGraph *filter) {
            std::stringstream ss(argv);
            std::string fileName;
            int libNum;
            int step;
            int mxVal;
            ss >> fileName >> libNum >> step >> mxVal;

            INFO("Weight statistic fileName=" << fileName << " libNum=" << libNum << " step=" << step << " mxVal="
                                              << mxVal);

            statistics::WeightStatistic cs;
            cs.calculateStatistic(filter, fileName, libNum, step, mxVal);
        }
    }
}
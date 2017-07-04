#include <Filter/Statistics/WeightDifStatistic.h>
#include "CommandWeightDifStatistic.h"

namespace filter {
    namespace commands {
        void CommandWeightDifStatistic::execute(std::string argv, State &state, Filter *filter) {
                std::stringstream ss(argv);
                std::string fileName;
                int libNum;
                int stepW, stepD;
                int mxW, mxD;
                ss >> fileName >> libNum >> stepW >> mxW >> stepD >> mxD;

                INFO("WeightDif statistic fileName=" << fileName << " libNum=" << libNum
                                                     << " stepW=" << stepW << " mxW=" << mxW
                                                     << " stepD=" << stepD << " mxD=" << mxD);

                statistics::WeightDifStatistic cs;
                cs.calculateStatistic(filter, fileName, libNum, stepW, mxW, stepD, mxD);
        }
    }
}
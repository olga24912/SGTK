#include <Filter/Statistics/TwoLibStatistic.h>
#include "CommandTwoLibStatistic.h"

namespace filter {
    namespace commands {

        void CommandTwoLibStatistic::execute(std::string argv, State &state, ContigGraph &graph) {
                std::stringstream ss(argv);
                std::string fileName;
                int libNum1;
                int step1;
                int mxVal1;
                int step1_2;
                int mxVal1_2;
                int libNum2;
                int step2;
                int mxVal2;
                ss >> fileName >> libNum1 >> step1 >> mxVal1 >> step1_2 >> mxVal1_2 >> libNum2 >> step2 >> mxVal2;

                INFO("two lib statistic fileName=" << fileName << " libNum1=" << libNum1 << " step1=" << step1
                                                   << " mxVal1="
                                                   << mxVal1
                                                   << " step1_2=" << step1_2 << " mxVal1_2=" << mxVal1_2 << " libNum2="
                                                   << libNum2
                                                   << " step2=" << step2 << " mxVal2=" << mxVal2);

                statistics::TwoLibStatistic cs;
                cs.calculateStatistic(&graph, fileName, libNum1, step1, mxVal1, step1_2, mxVal1_2, libNum2, step2,
                                      mxVal2);
        }
    }
}
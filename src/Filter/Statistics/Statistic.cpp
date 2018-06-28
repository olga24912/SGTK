#include <iostream>
#include "Statistic.h"

namespace filter {
    namespace statistics {
        void Statistic::printStatistic(char *str, int *cnt) {
            int sum = 0;
            for (int i = 0; i < 7; ++i) {
                sum += cnt[i];
            }
            if (sum < 1) return;
            if (sum == 0) sum = 1;
            double TP = cnt[0] * 1.0 / sum, FP = (cnt[3] + cnt[4] + cnt[5]) * 1.0 / sum;

            printf("%s\tTP - %lf; \t FP - %lf \t "
                           "OK - %d; \t OVERLAP - %d; \t PART_ALIG - %d;\t"
                           " BIG_DIST - %d; \t WRONG_ORDER - %d;\t  DIF_CHR - %d; \t NA - %d;\n",
                   str, TP, FP,
                   cnt[0], cnt[1], cnt[2], cnt[3],
                   cnt[4], cnt[5], cnt[6]);
        }
    }
}

#ifndef SCAFFOLDER_STATISTIC_H
#define SCAFFOLDER_STATISTIC_H


#include "InfoAboutContigsAlig.h"

class Statistic {
protected:
    enum ErrorType {OK, OVERLAP, PART_ALIG, BIG_DIST, WRONG_ORDER, DIF_CHR, NA};
    const int MAX_DIST = 1000000;
    const int MIN_OVERLAP = 100;
    ErrorType isCorrectEdge(InfoAboutContigsAlig& aligInfo, Filter * filter, int e);
};


#endif //SCAFFOLDER_STATISTIC_H

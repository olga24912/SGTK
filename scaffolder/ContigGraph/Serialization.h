//
// Created by olga on 23.10.16.
//

#ifndef SCAFFOLDER_SERIALIZATION_H
#define SCAFFOLDER_SERIALIZATION_H

#include "ContigGraph.h"

using namespace std;

class Serialization {
public:
    static void write(ContigGraph * g, string fileName);
    static ContigGraph read(string fileName);
};


#endif //SCAFFOLDER_SERIALIZATION_H

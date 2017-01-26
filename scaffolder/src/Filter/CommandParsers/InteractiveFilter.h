/*//
// Created by olga on 24.10.16.
//

#ifndef SCAFFOLDER_INTERACTIVEFILTER_H
#define SCAFFOLDER_INTERACTIVEFILTER_H

#include "ContigGraph/ContigGraph.h"
#include "ContigGraph/Serialization.h"
#include "Filter/Scaffolder/ScafSimplePath.h"
#include "Filter/Tools/SystemTools.h"
#include "Filter/Tools/Utils.h"
#include "ContigGraphPrinter.h"

class InteractiveFilter {
private:

    ContigGraph g;

    void readConfig();
    bool handlingRequest(istream &in);
    void handlingWriteRequest(string req, istream &in);

public:
};


#endif //SCAFFOLDER_INTERACTIVEFILTER_H
*/
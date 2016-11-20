//
// Created by olga on 19.11.16.
//

#ifndef SCAFFOLDER_REFERENCEGRAPHBUILDER_H
#define SCAFFOLDER_REFERENCEGRAPHBUILDER_H


#include "GraphBuilder.h"
#include <seqan/seq_io.h>
#include "../Tools/WorkWithOtherTools.h"
#include "../Tools/SeqanUtils.h"


class ReferenceGraphBuilder: public GraphBuilder {
private:
    const int MIN_CONTIG = 500;
    string refContigFileName, queryContigsFileName;

    map <string, int> contigsId;
    vector<string> contigsName;

    void generateVertex();
    void createGraph(string fileName);
public:
    void evaluate();

    void setRefFileName(const string &refContigFileName) {
        ReferenceGraphBuilder::refContigFileName = refContigFileName;
    }

    void setQueryFileName(const string &queryContigsFileName) {
        ReferenceGraphBuilder::queryContigsFileName = queryContigsFileName;
    }
};


#endif //SCAFFOLDER_REFERENCEGRAPHBUILDER_H

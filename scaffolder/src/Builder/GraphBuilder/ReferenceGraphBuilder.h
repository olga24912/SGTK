//
// Created by olga on 19.11.16.
//

#ifndef SCAFFOLDER_REFERENCEGRAPHBUILDER_H
#define SCAFFOLDER_REFERENCEGRAPHBUILDER_H


#include "GraphBuilder.h"
#include <seqan/seq_io.h>
#include "Builder/Tools/SystemAlignmentTools.h"
#include "Builder/Tools/SeqanUtils.h"


class ReferenceGraphBuilder: public GraphBuilder {
private:
    const int MIN_CONTIG = 500;
    string refContigFileName, queryContigsFileName;
    string tsvFileName;

    struct alignmentInfo {
        int sr;
        int er;
        int sq;
        int eq;
        string contigName;

        alignmentInfo(){}
        alignmentInfo(int sr, int er, int sq, int eq, string name): sr(sr), er(er), sq(sq), eq(eq), contigName(name) {}

        bool operator < (alignmentInfo b) {
            return (sr < b.sr);
        }
    };

    map <string, int> contigsId;
    vector<string> contigsName;

    void generateVertex();
    void createGraph(map<string, vector<alignmentInfo>> contigsAlignment);
    map<string, vector<alignmentInfo>> parseCoordFile(string fileName);
    map<string, vector<alignmentInfo>> parseTSVFile(string fileName);
    virtual string getLibColor();
public:
    void evaluate();

    void setRefFileName(const string &refContigFileName) {
        ReferenceGraphBuilder::refContigFileName = refContigFileName;
    }

    void setQueryFileName(const string &queryContigsFileName) {
        ReferenceGraphBuilder::queryContigsFileName = queryContigsFileName;
    }

    void setTsvFileName(const string &tsvFileName) {
        ReferenceGraphBuilder::tsvFileName = tsvFileName;
    }
};


#endif //SCAFFOLDER_REFERENCEGRAPHBUILDER_H

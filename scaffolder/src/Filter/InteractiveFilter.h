//
// Created by olga on 24.10.16.
//

#ifndef SCAFFOLDER_INTERACTIVEFILTER_H
#define SCAFFOLDER_INTERACTIVEFILTER_H

#include "ContigGraph/ContigGraph.h"
#include "ContigGraph/Serialization.h"
#include "Filter/Scaffolder/ScafSimplePath.h"
#include "Tools/SystemTools.h"

class InteractiveFilter {
private:
    static const string UPLOAD_GRAPH;
    static const string MIN_EDGE_WEIGHT;
    static const string MIN_CONTIG_LEN;
    static const string WRITE_FULL;
    static const string WRITE_LOCAL;
    static const string WRITE_LOCAL_VERT_IN_SEG;
    static const string WRITE_ALL_LOCAL;
    static const string WRITE_BIG_COMP;
    static const string WRITE_SPLIT_BIG_COMP;
    static const string MERGE_SIMPLE_PATH;
    static const string WRITE_LOCAL_ALONG_PATH;
    static const string SET_IGNORE;
    static const string RESET_IGNORE;
    static const string EXIT;

    static const string CONFIG_FILE;

    ContigGraph g;

    void readConfig();
    bool handlingRequest(istream &in);
    void handlingWriteRequest(string req, istream &in);

    struct State {
        enum StateName {NEW, LOCAL};
        StateName name = StateName::NEW;
        string fileName = "";
        int dist = 0;
    };

    State state;
public:
    void main();
};


#endif //SCAFFOLDER_INTERACTIVEFILTER_H
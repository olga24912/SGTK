#ifndef SCAFFOLDER_SCAFSIMPLEPATH_H
#define SCAFFOLDER_SCAFSIMPLEPATH_H


#include "ContigGraph/ContigGraph.h"
#include <seqan/seq_io.h>
#include "Tools/SeqanUtils.h"

using namespace std;
using namespace seqan;

class ScafSimplePath {
private:
    const int GAP_SIZE = 100;

    ContigGraph* graph;
    string contigFile;

    vector<string> contigs;
    vector<int> next;

    void saveContigs();
    void findPaths();
    void writeNewContigs(string out);

    string createRevCompl(string s);
    void dfsPath(int v, vector<int> &next, vector<int> &used);
public:
    void evaluate(ContigGraph* graph, string contigFile, string out);
};


#endif //SCAFFOLDER_SCAFSIMPLEPATH_H

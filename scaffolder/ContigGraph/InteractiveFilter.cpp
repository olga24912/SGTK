//
// Created by olga on 24.10.16.
//

#include "InteractiveFilter.h"
#include "ContigGraphPrinter.h"



const string InteractiveFilter::UPLOAD_GRAPH = "uploadGraph";
const string InteractiveFilter::MIN_EDGE_WEIGHT = "minEdgeW";
const string InteractiveFilter::MIN_CONTIG_LEN = "minContig";
const string InteractiveFilter::WRITE_FULL = "writeFull";
const string InteractiveFilter::WRITE_LOCAL = "writeLocal";
const string InteractiveFilter::EXIT = "exit";

/*
 * commands:
 * uploadGraph <filename>
 * minEdgeW <libNum> <weight>
 * minContig <len>
 * writeFull <fileName>
 * writeLocal <vertexID> <dist> <fileName>
 * exit
 */
void InteractiveFilter::main() {
    ContigGraph g;
    while(true) {
        string s;
        cin >> s;
        if (s == UPLOAD_GRAPH) {
            string fileName;
            cin >> fileName;
            g = Serialization::read(fileName);
        } else if (s == MIN_EDGE_WEIGHT) {
            int libn;
            int w;
            cin >> libn >> w;
            g.setMinEdgeWeightForLib(libn, w);
        } else if (s == MIN_CONTIG_LEN) {
            int len;
            cin >> len;
            g.filterByContigLen(len);
        } else if (s == WRITE_FULL) {
            string fileName;
            cin >> fileName;
            ContigGraphPrinter::writeFullGraphDotFormat(&g, fileName);
        } else if (s == WRITE_LOCAL) {
            int v;
            int dist;
            string fileName;
            cin >> v >> dist >> fileName;
            ContigGraphPrinter::writeLocalGraph(&g, dist, v, fileName);
        } else if (s == EXIT) {
            return;
        }
    }
}

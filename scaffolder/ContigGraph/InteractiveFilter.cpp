#include "InteractiveFilter.h"
#include "ContigGraphPrinter.h"

const string InteractiveFilter::UPLOAD_GRAPH = "uploadGraph";
const string InteractiveFilter::MIN_EDGE_WEIGHT = "minEdgeW";
const string InteractiveFilter::MIN_CONTIG_LEN = "minContig";
const string InteractiveFilter::WRITE_FULL = "writeFull";
const string InteractiveFilter::WRITE_LOCAL = "writeLocal";
const string InteractiveFilter::WRITE_ALL_LOCAL = "writeAllLocal";
const string InteractiveFilter::WRITE_LOCAL_VERT_IN_SEG = "writeLocalSeg";
const string InteractiveFilter::WRITE_BIG_COMP = "writeBig";
const string InteractiveFilter::WRITE_SPLIT_BIG_COMP = "writeSB";
const string InteractiveFilter::MERGE_SIMPLE_PATH = "mergeSimplePath";
const string InteractiveFilter::EXIT = "exit";

/*
 * commands:
 * uploadGraph <filename>
 * minEdgeW <libNum> <weight>
 * minContig <len>
 * writeFull <fileName>
 * writeLocal <vertexID> <dist> <fileName>
 * writeAllLocal <dist>
 * writeLocalSeg <vertexIDStart> <vertexIDFinish> <dist> <fileName>
 * writeBig <size> <fileName>
 * writeSB <size> <prefixFileName>
 * mergeSimplePath <contigsFileName> <outFileName>
 * exit
 */

void InteractiveFilter::main() {
    ContigGraph g;
    while(true) {
        cout << "ok" << endl;
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
        } else if (s == WRITE_ALL_LOCAL) {
            int dist;
            cin >> dist;
            ContigGraphPrinter::writeAllLocalGraphDotFormat(&g, dist);
        } else if (s == WRITE_LOCAL_VERT_IN_SEG) {
            int vb, ve;
            int dist;
            string fileName;
            cin >> vb >> ve >> dist >> fileName;
            ContigGraphPrinter::writeLocalSegGraph(&g, dist, vb, ve, fileName);
        } else if (s == WRITE_BIG_COMP) {
            int size;
            string fileName;
            cin >> size >> fileName;
            ContigGraphPrinter::writeBigComponent(&g, size, fileName);
        } else if (s == WRITE_SPLIT_BIG_COMP) {
            int size;
            string fileName;
            cin >> size >> fileName;
            ContigGraphPrinter::writeSplitBigComponent(&g, size, fileName);
        } else if (s == MERGE_SIMPLE_PATH) {
            string fileNameIn;
            string fileNameOut;
            cin >> fileNameIn >> fileNameOut;
            ScafSimplePath ssp;
            ssp.evaluate(&g, fileNameIn, fileNameOut);
        } else if (s == EXIT) {
            return;
        }
    }
}

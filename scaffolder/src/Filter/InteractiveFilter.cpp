#include "InteractiveFilter.h"

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
const string InteractiveFilter::WRITE_LOCAL_ALONG_PATH = "writeAlongPath";
const string InteractiveFilter::SET_IGNORE = "setIgnore";
const string InteractiveFilter::RESET_IGNORE = "resetIgnore";
const string InteractiveFilter::EXIT = "exit";

const string InteractiveFilter::CONFIG_FILE = "filter_config";

/*
 * commands:
 * uploadGraph <filename>
 * minEdgeW <libNum> <weight>
 * minContig <len>
 * mergeSimplePath <contigsFileName> <outFileName>
 * setIgnore <vertexIdStart> <vertexIdFinish>
 * resetIgnore
 * exit
 * writeFull <fileName>
 * writeLocal <fileName> <vertexID> <dist>
 * writeAllLocal <fileName> <dist>
 * writeLocalSeg <fileName> <vertexIDStart> <vertexIDFinish> <dist>
 * writeBig <fileName> <size>
 * writeSB <prefixFileName> <size>
 * writeAlongPath <prefixFileName> <libNum> <dist> <minRefPathSize>
 */

void InteractiveFilter::main() {
    readConfig();

    while(handlingRequest(cin)) {
        cout << "ok" << endl;
    }
}

void InteractiveFilter::readConfig() {
    ifstream in(CONFIG_FILE);

    while (handlingRequest(in));
}

bool InteractiveFilter::handlingRequest(istream &in) {
    string s;
    if (!(in >> s)) {
        return false;
    }
    if (s == UPLOAD_GRAPH) {
        string fileName;
        in >> fileName;
        g = Serialization::read(fileName);
    } else if (s == MIN_EDGE_WEIGHT) {
        int libn;
        int w;
        in >> libn >> w;
        g.setMinEdgeWeightForLib(libn, w);
    } else if (s == MIN_CONTIG_LEN) {
        int len;
        in >> len;
        g.filterByContigLen(len);
    } else if (s == MERGE_SIMPLE_PATH) {
        string fileNameIn;
        string fileNameOut;
        in >> fileNameIn >> fileNameOut;
        ScafSimplePath ssp;
        ssp.evaluate(&g, fileNameIn, fileNameOut);
    } else if (s == SET_IGNORE) {
        int vs;
        int vf;
        in >> vs >> vf;
        for (int i = vs; i <= vf; ++i) {
            g.setIgnore(i);
        }
    } else if (s == RESET_IGNORE) {
        for (int i = 0; i < g.getVertexCount(); ++i) {
            if (g.isIgnore(i)) {
                g.setIgnore(i);
            }
        }
    } else if (s == EXIT) {
            return false;
    } else {
        handlingWriteRequest(s, in);
    }
    return true;
}

void InteractiveFilter::handlingWriteRequest(string s, istream& in) {
    string fileName;
    if (!Utils::isInt(s)) {
        in >> fileName;
    }
    if (s == WRITE_FULL) {
        ContigGraphPrinter::writeFullGraphDotFormat(&g, fileName);
    } else if (s == WRITE_LOCAL) {
        int v;
        int dist;
        in >> v >> dist;
        ContigGraphPrinter::writeLocalGraph(&g, dist, v, fileName);
        state.name = State::LOCAL;
        state.fileName = fileName;
        state.dist = dist;
    } else if (s == WRITE_ALL_LOCAL) {
        int dist;
        in >> dist;
        ContigGraphPrinter::writeAllLocalGraphDotFormat(&g, dist);
    } else if (s == WRITE_LOCAL_VERT_IN_SEG) {
        int vb, ve;
        int dist;
        in >> vb >> ve >> dist;
        ContigGraphPrinter::writeLocalSegGraph(&g, dist, vb, ve, fileName);
    } else if (s == WRITE_BIG_COMP) {
        int size;
        cin >> size;
        ContigGraphPrinter::writeBigComponent(&g, size, fileName);
    } else if (s == WRITE_SPLIT_BIG_COMP) {
        int size;
        cin >> size;
        ContigGraphPrinter::writeSplitBigComponent(&g, size, fileName);
    } else if (s == WRITE_LOCAL_ALONG_PATH) {
        int libId;
        int dist;
        int minSize;
        in >> libId >> dist >> minSize;
        ContigGraphPrinter::writeAlongPath(&g, libId, dist, minSize, fileName);
    } else if (Utils::isInt(s)) {
        if (state.name == State::LOCAL) {
            fileName = state.fileName;
            int dist = state.dist;
            int v = Utils::strToInt(s);
            ContigGraphPrinter::writeLocalGraph(&g, dist, v, fileName + "_" + s);
        }
    }

    SystemTools::showDotFile(fileName);
}

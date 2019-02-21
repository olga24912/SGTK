#include "GraphBuilder.h"

namespace builder {
    namespace graph_builder {
        std::string GraphBuilder::dna5ToString(seqan::Dna5 *seq, int len) {
            std::string res;

            std::string intToDNAChar = "ACGTN";

            for (int i = 0; i < len; ++i) {
                res += intToDNAChar[ordValue(seq[i])];
            }
            return res;
        }

        std::string GraphBuilder::cutReadName(seqan::BamAlignmentRecord read) {
            std::string readName = std::string(seqan::toCString(read.qName));
            if (readName.size() > 1) {
                if ((readName[readName.size() - 2] == '/' ||
                     readName[readName.size() - 2] == '_') &&
                    (readName[readName.size() - 1] == '2' ||
                     readName[readName.size() - 1] == '1')) {
                    readName.resize(readName.size() - 2);
                }
            }
            return readName;
        }

        void GraphBuilder::setGraph(ContigGraph *graph) {
            DEBUG("setGraph");
            GraphBuilder::graph = graph;
            graph->newLib(libName, getLibColor(), getLibType());
        }

        void GraphBuilder::setLibName(std::string libName) {
            DEBUG("setLibName libName=" << libName);
            this->libName = libName;
        }

        void GraphBuilder::setContigFile(const std::string& file_name) {
            contigFileName = file_name;
        }
    }
}

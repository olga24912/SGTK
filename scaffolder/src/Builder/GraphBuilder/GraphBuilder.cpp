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

        std::string GraphBuilder::colorToString(int *color) {
            std::string res = "#";
            for (int i = 0; i < 3; ++i) {
                if (color[i] / 16 < 10) {
                    res += (color[i] / 16) + '0';
                } else {
                    res += (color[i] / 16) - 10 + 'a';
                }

                if (color[i] % 16 < 10) {
                    res += (color[i] % 16) + '0';
                } else {
                    res += (color[i] % 16) - 10 + 'a';
                }
            }
            TRACE("colorToString color=(" << color[0] << " " << color[1] << " " << color[2] << ") : " << res);
            return res;
        }

        void GraphBuilder::setSamFileWriter() {
            using namespace sam_file_writer;
            DEBUG("setSamFileWriter");
            this->samFileWriter = SamFileWriteEdge(path);
        }
    }
}

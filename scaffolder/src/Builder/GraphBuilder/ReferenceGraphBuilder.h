#ifndef SCAFFOLDER_REFERENCEGRAPHBUILDER_H
#define SCAFFOLDER_REFERENCEGRAPHBUILDER_H

#include <seqan/seq_io.h>
#include "GraphBuilder.h"
#include "Builder/Tools/SystemAlignmentTools.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
//It understand contig order by alignment contigs on reference assembly
//It is nessery put ((ref and query files name) or tsv file name)
        class ReferenceGraphBuilder : public GraphBuilder {
        private:
            int minContigLen = 500;

            std::string refContigFileName;
            std::string queryContigsFileName;
            std::string tsvFileName;

            struct alignmentInfo {
                int sr;
                int er;
                int sq;
                int eq;
                std::string contigName;

                alignmentInfo() {}

                alignmentInfo(int sr, int er, int sq, int eq, std::string name) : sr(sr), er(er), sq(sq), eq(eq),
                                                                                  contigName(name) {}

                bool operator<(alignmentInfo b) {
                    return (sr < b.sr);
                }
            };

            std::map<std::string, int> contigsId;
            std::vector<std::string> contigsName;

            void generateVertex();

            void createGraph(std::map<std::string, std::vector<alignmentInfo>> contigsAlignment);

            std::map<std::string, std::vector<alignmentInfo>> parseCoordFile(std::string fileName);

            std::map<std::string, std::vector<alignmentInfo>> parseTSVFile(std::string fileName);

            virtual std::string getLibColor();

        protected:
            ContigGraph::Lib::Type getLibType() override;

        public:
            void evaluate();

            //file with reference assembly
            void setRefFileName(const std::string &refContigFileName) {
                ReferenceGraphBuilder::refContigFileName = refContigFileName;
            }

            //file with contigs some assembly.
            void setQueryFileName(const std::string &queryContigsFileName) {
                ReferenceGraphBuilder::queryContigsFileName = queryContigsFileName;
            }

            //file with information about alignment in tsv format
            void setTsvFileName(const std::string &tsvFileName) {
                ReferenceGraphBuilder::tsvFileName = tsvFileName;
            }

            // set barrier contig len. Contigs with smaller len will be ignore.
            void setMinContigLen(int minContigLen);

        private:
            DECL_LOGGER("ReferenceGraphBuilder");
        };
    }
}

#endif //SCAFFOLDER_REFERENCEGRAPHBUILDER_H

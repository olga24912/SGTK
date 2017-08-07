#include "RNASplitReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        void RNASplitReadGraphBuilder::evaluate() {
            using namespace tools;
            using namespace reads_splitter;
            INFO("start RNASplitReadGraph build");
            SystemAlignmentTools wwot;
            ReadsSplitter50 rs;
            SplitterByUnmappedEnd su;

            wwot.alignmentRNA(refFileName, rnaReadsFileName, "rna.sam", path);

            std::string unmappedName;
            if (rnaReadsFileName[rnaReadsFileName.size() - 1] == 'q') {
                unmappedName = path + "/Unmapped.fastq";
            } else {
                unmappedName = path + "/Unmapped.fasta";
            }
            std::string command = "mv " + path + "/Unmapped.out.mate1 " + unmappedName;
            system(command.c_str());

            rs.splitReads(unmappedName, path + "/cutPartReads1.fasta", path + "/cutPartReads2.fasta");
            su.splitReads(path + "/rna.sam", path + "/short1.fasta", path + "/short2.fasta");

            handlingPairReads(path + "/cutPartReads1.fasta", path + "/cutPartReads2.fasta", libName + "-50-50");
            graph->setLib(ContigGraph::Lib::RNA_SPLIT_50);
            handlingPairReads(path + "/short1.fasta", path + "/short2.fasta", libName + "-long-short");
            graph->setLib(ContigGraph::Lib::RNA_SPLIT_30);
            INFO("finish RNASplitReadGraph build");
        }

        void RNASplitReadGraphBuilder::setRefFileName(std::string refFileName) {
            TRACE("setRefFileName");
            RNASplitReadGraphBuilder::refFileName = refFileName;
        }

        void RNASplitReadGraphBuilder::setRnaReadFileName(std::string rnaReadsFileName) {
            TRACE("setRnaReadFileName");
            RNASplitReadGraphBuilder::rnaReadsFileName = rnaReadsFileName;
        }

        void RNASplitReadGraphBuilder::handlingPairReads(std::string file1, std::string file2, std::string libN) {
            INFO("start handle pair reads file1=" << file1 << " file2=" << file2 << " libName=" << libN);
            tools::SystemAlignmentTools wwot;

            wwot.alignmentRNA(refFileName, file1, "rna1.sam", path);
            wwot.alignmentRNA(refFileName, file2, "rna2.sam", path);

            RNAPairReadGraphBuilder gb;

            gb.setLibName(libN, path);
            gb.setFileName1(path + "/rna1.sam");
            gb.setFileName2(path + "/rna2.sam");
            gb.setOneSideReadFlag(true);
            gb.setGraph(graph);

            gb.evaluate();

            INFO("finish handle pair reads");
        }

        std::string RNASplitReadGraphBuilder::getLibColor() {
            TRACE("getLibColor");
            int cntRB = 100 + rand() % 150;
            int color[3] = {cntRB, rand() % (cntRB / 2), cntRB};
            return colorToString(color);
        }

        void RNASplitReadGraphBuilder::setGraph(ContigGraph *graph) {
            TRACE("setGraph");
            GraphBuilder::graph = graph;
        }

        void RNASplitReadGraphBuilder::setSamFileWriter() {}

        ContigGraph::Lib::Type RNASplitReadGraphBuilder::getLibType() {
            TRACE("getLibType");
            return ContigGraph::Lib::RNA_SPLIT_50;
        }
    }
}
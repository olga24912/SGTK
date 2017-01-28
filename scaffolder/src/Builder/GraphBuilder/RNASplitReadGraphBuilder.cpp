#include "RNASplitReadGraphBuilder.h"

void RNASplitReadGraphBuilder::evaluate() {
    SystemAlignmentTools wwot;
    ReadsSplitter50 rs;
    SplitterByUnmappedEnd su;

    wwot.alignmentRNA(refFileName, rnaReadsFileName, "rna.sam", path);
    string unmappedName;
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
    handlingPairReads(path + "/short1.fasta", path + "/short2.fasta", libName + "-long-short");
}

void RNASplitReadGraphBuilder::setRefFileName(string refFileName) {
    RNASplitReadGraphBuilder::refFileName = refFileName;
}

void RNASplitReadGraphBuilder::setRnaReadFileName(string rnaReadsFileName) {
    RNASplitReadGraphBuilder::rnaReadsFileName = rnaReadsFileName;
}

void RNASplitReadGraphBuilder::handlingPairReads(string file1, string file2, string libN) {
    SystemAlignmentTools wwot;

    wwot.alignmentRNA(refFileName, file1, "rna1.sam", path);
    wwot.alignmentRNA(refFileName, file2, "rna2.sam", path);

    RNAPairReadGraphBuilder gb;

    gb.setLibName(libN, path);
    gb.setFileName1(path + "/rna1.sam");
    gb.setFileName2(path + "/rna2.sam");
    gb.setOneSideReadFlag(true);
    gb.setGraph(graph);

    gb.evaluate();
}

string RNASplitReadGraphBuilder::getLibColor() {
    int cntRB = 100 + rand()%150;
    int color[3] = {cntRB, rand()%(cntRB/2), cntRB};
    return colorToString(color);
}

void RNASplitReadGraphBuilder::setGraph(ContigGraph *graph) {
    GraphBuilder::graph = graph;
}

void RNASplitReadGraphBuilder::setSamFileWriter() {}

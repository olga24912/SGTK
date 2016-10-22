//
// Created by olga on 08.10.16.
//

#include "RNASplitReadGraphBuilder.h"
#include "../Tools/WorkWithOtherTools.h"
#include "RNAPairReadGraphBuilder.h"
#include "../ReadsSplitter/ReadsSplitter.h"
#include "../ReadsSplitter/ReadsSplitter50.h"
#include "../ReadsSplitter/SplitterByUnmappedEnd.h"

void RNASplitReadGraphBuilder::evaluate() {
    WorkWithOtherTools wwot;
    ReadsSplitter50 rs;
    SplitterByUnmappedEnd su;

    wwot.alignmentRNA(refFileName, rnaReadsFileName, "rna.sam");
    system("mv Unmapped.out.mate1 Unmapped.fasta");
    rs.splitReads("Unmapped.fasta", "cutPartReads1.fasta", "cutPartReads2.fasta");
    su.splitReads("rna.sam", "short1.fasta", "short2.fasta");

    handlingPairReads("cutPartReads1.fasta", "cutPartReads2.fasta");
    handlingPairReads("short1.fasta", "short2.fasta");

    graph->filterByContigLen(minContigLen);
}

void RNASplitReadGraphBuilder::setRefFileName(string refFileName) {
    RNASplitReadGraphBuilder::refFileName = refFileName;
}

void RNASplitReadGraphBuilder::setRnaReadFileName(string rnaReadsFileName) {
    RNASplitReadGraphBuilder::rnaReadsFileName = rnaReadsFileName;
}

void RNASplitReadGraphBuilder::handlingPairReads(string file1, string file2) {
    WorkWithOtherTools wwot;

    wwot.alignmentRNA(refFileName, file1, "rna1.sam");
    wwot.alignmentRNA(refFileName, file2, "rna2.sam");

    RNAPairReadGraphBuilder gb;

    gb.setFileName1("rna1.sam");
    gb.setFileName2("rna2.sam");
    gb.setOneSideReadFlag(true);
    gb.setGraph(graph);
    gb.setMinContigLen(minContigLen);
    gb.setMinEdgeWight(minEdgeWight);

    gb.evaluate();
    graph->filterByEdgeWeight(minEdgeWight);
}

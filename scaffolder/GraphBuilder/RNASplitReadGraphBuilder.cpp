//
// Created by olga on 08.10.16.
//

#include "RNASplitReadGraphBuilder.h"
#include "../Tools/WorkWithOtherTools.h"
#include "../Tools/FastaToolsIn.h"
#include "../ReadsAlignment.h"
#include "RNAPairReadGraphBuilder.h"
#include "../Tools/ReadsSplitter.h"

void RNASplitReadGraphBuilder::evaluate() {
    WorkWithOtherTools wwot;
    ReadsSplitter rs;
    ReadsAlignment ra;


    wwot.alignmentRNA(refFileName, rnaReadsFileName, "rna.sam");
    rs.findAndSplitNotAlignmentReads(rnaReadsFileName, "rna.sam", "cur2reads.fasta");
    wwot.alignmentRNA(refFileName, "cur2reads.fasta", "rnaSplit.sam");
    /*ra.alignment(refFileName, "rnaSplit.sam", "cutPartReads1.fasta", "cutPartReads2.fasta");


    wwot.alignmentRNA(refFileName, "cutPartReads1.fasta", "rna1.sam");
    wwot.alignmentRNA(refFileName, "cutPartReads2.fasta", "rna2.sam");

    RNAPairReadGraphBuilder gb;

    gb.setFileName1("rna1.sam");
    gb.setFileName2("rna2.sam");
    gb.setGraph(graph);
    gb.setMinContigLen(minContigLen);
    gb.setMinEdgeWight(minEdgeWight);

    gb.evaluate();*/
}

void RNASplitReadGraphBuilder::setRefFileName(string refFileName) {
    RNASplitReadGraphBuilder::refFileName = refFileName;
}

void RNASplitReadGraphBuilder::setRnaReadFileName(string rnaReadsFileName) {
    RNASplitReadGraphBuilder::rnaReadsFileName = rnaReadsFileName;
}

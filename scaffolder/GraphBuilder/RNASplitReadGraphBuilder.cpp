//
// Created by olga on 08.10.16.
//

#include "RNASplitReadGraphBuilder.h"
#include "../Tools/WorkWithOtherTools.h"

void RNASplitReadGraphBuilder::evaluate() {
    WorkWithOtherTools wwot;

    wwot.alignmentRNA(refFileName, rnaReadsFileName, "rna.sam");
}

void RNASplitReadGraphBuilder::setRefFileName(string refFileName) {
    RNASplitReadGraphBuilder::refFileName = refFileName;
}

void RNASplitReadGraphBuilder::setRnaReadFileName(string rnaReadsFileName) {
    RNASplitReadGraphBuilder::rnaReadsFileName = rnaReadsFileName;
}

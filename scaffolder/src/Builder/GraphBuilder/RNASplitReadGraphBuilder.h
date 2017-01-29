//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H
#define SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include "Builder/Tools/SystemAlignmentTools.h"
#include "RNAPairReadGraphBuilder.h"
#include "Builder/ReadsSplitter/ReadsSplitter.h"
#include "Builder/ReadsSplitter/ReadsSplitter50.h"
#include "Builder/ReadsSplitter/SplitterByUnmappedEnd.h"


// find connection between contig by read who split on exons on
// two deferent contigs.
class RNASplitReadGraphBuilder : public GraphBuilder {
private:
    std::string refFileName;
    std::string rnaReadsFileName;

    virtual std::string getLibColor();

protected:
    void setSamFileWriter() override;
    void handlingPairReads(std::string file1, std::string file2, std::string linN);

public:
    void evaluate();
    // set name of fie with ref genom in fasta format.
    void setRefFileName(std::string refFileName);

    // set name of  file with rna reads in fasta format.
    void setRnaReadFileName(std::string rnaReadsFileName);

    void setGraph(ContigGraph *graph);
};


#endif //SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

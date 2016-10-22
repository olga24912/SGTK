//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H
#define SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

#include "GraphBuilder.h"

/*
 * find connection between contig by read who split on exons on
 * two deferent contigs.
 */
class RNASplitReadGraphBuilder : public GraphBuilder {
private:
    string refFileName;
    string rnaReadsFileName;

public:
    void evaluate();
    /*
     * set name of fie with ref genom in fasta format.
     */
    void setRefFileName(string refFileName);

    /*
     * set name of  file with rna reads in fasta format.
     */
    void setRnaReadFileName(string rnaReadsFileName);

    void handlingPairReads(string file1, string file2, string linN);
};


#endif //SCAFFOLDER_RNASPLITREADGRAPHBUILDER_H

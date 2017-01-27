//
// Created by olga on 08.10.16.
//

#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H


#include "ContigGraph/ContigGraph.h"
#include "Builder/SamFileWriter/SamFileWriteEdge.h"

/*
 * main class for generate conection between contigs.
 */
class GraphBuilder {
protected:
    ContigGraph* graph;
    SamFileWriteEdge samFileWriter;

    string libName;
    int minContigLen = 0;

    virtual string getLibColor() = 0;
    static string colorToString(int color[3]);
public:
    //fun that need to call for add conection between contigs;
    virtual void evaluate() = 0;

    // set barrier contig len. Contigs with smaller len will be ignore.
    void setMinContigLen(int minContigLen);

    virtual void setGraph(ContigGraph* graph);

    void setLibName(string libName);
    void setSamFileWriter(SamFileWriteEdge fileWriter);
};


#endif //SCAFFOLDER_GRAPHBUILDER_H

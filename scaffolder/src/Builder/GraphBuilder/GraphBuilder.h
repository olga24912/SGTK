#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H

#include <Builder/SamFileWriter/SamFileWriterEdgeFictive.h>
#include "ContigGraph/ContigGraph.h"
#include "Builder/SamFileWriter/SamFileWriteEdge.h"

// main class for generate conection between contigs.
class GraphBuilder {
protected:
    std::string path; // path/to/the/dir where this class will generate files.

    ContigGraph* graph; //generated graph
    SamFileWriteEdge samFileWriter = SamFileWriterEdgeFictive(); //for write info about edge's reads

    std::string libName; //name of the lib

    virtual std::string getLibColor() = 0; //return color for current lib
    virtual ContigGraph::Lib::Type getLibType() = 0;
    static std::string colorToString(int color[3]); // translate color like array {255, 0, 255} to  string "#ff00ff"
public:
    //fun that need to call for add conection between contigs;
    virtual void evaluate() = 0;

    //set graph, and change lib in it to new.
    // need to be coll after setLibName
    virtual void setGraph(ContigGraph* graph);

    //set libName
    //create new dir: path/libName
    void setLibName(std::string libName, std::string path);

    virtual void setSamFileWriter(); // init samFileWriter
};


#endif //SCAFFOLDER_GRAPHBUILDER_H

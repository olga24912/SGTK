#ifndef SCAFFOLDER_SIMPLEWRITER_H
#define SCAFFOLDER_SIMPLEWRITER_H

#include <fstream>
#include <algorithm>
#include <Filter/Writers/FileValidator/FileValidator.h>
#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include <Filter/AligInfo/InfoAboutContigsAlig.h>
#include <Filter/ContigGraph/ContigGraph.h>
#include "Filter/Writers/GraphSplitter.h"

namespace filter {
    namespace writers {
        using namespace contig_graph;
//class for write graph in dot format.
        class DotWriter {
        protected:
            ContigGraph *graph; //graph for writing
            GraphSplitter graphSplitter; //split graph on small parts
            FileValidator *validator = new ValidatorNotPathWithAllLib();
            std::string coordFile = "";
            alig_info::InfoAboutContigsAlig aligInfo;

            virtual void writeOneVertex(int v, int isColored, std::ofstream &out);

            virtual void writeOneEdge(int e, std::ofstream &out);

            virtual void writeOneVertexSet(std::vector<int> vert, std::string fileName);

        public:
            DotWriter() {}

            DotWriter(ContigGraph *filter) : graph(filter) {}

            DotWriter(ContigGraph *filter, FileValidator *validator, int maxVert, int maxEdge) :
                    graph(filter), validator(validator) {
                graphSplitter = GraphSplitter(maxVert, maxEdge);
            }


            DotWriter(ContigGraph *filter, FileValidator *validator, int maxVert, int maxEdge, std::string coordFile) :
                    graph(filter), validator(validator), coordFile(coordFile) {
                graphSplitter = GraphSplitter(maxVert, maxEdge);
                if (coordFile != "") {
                    aligInfo.parseCoordFile(filter, coordFile);
                }
            }

            //write this set of vertex in files with prefix fileName
            virtual void writeVertexSet(std::vector<int> vert, std::string fileName);

        protected:
            DECL_LOGGER("DotWriter");

            void findVertWithOtherEdges(const std::vector<int> &vert, std::vector<bool> &hasOtherEdge,
                                std::vector<std::pair<int, int>> &weightEdge) const;
        };
    }
}

#endif //SCAFFOLDER_SIMPLEWRITER_H

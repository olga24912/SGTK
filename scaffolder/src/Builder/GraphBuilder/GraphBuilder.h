#ifndef SCAFFOLDER_GRAPHBUILDER_H
#define SCAFFOLDER_GRAPHBUILDER_H

#include <seqan/bam_io.h>
#include "Builder/ContigGraph/ContigGraph.h"

namespace builder {
    namespace graph_builder {
        using namespace contig_graph;
// main class for generate conection between contigs.
        class GraphBuilder {
        protected:
            std::string path; // path/to/the/dir where this class will generate files.

            ContigGraph *graph; //generated graph

            std::string libName; //name of the lib

            virtual std::string getLibColor() = 0; //return color for current lib
            virtual ContigGraph::Lib::Type getLibType() = 0;

            static std::string
            colorToString(int color[3]); // translate color like array {255, 0, 255} to  string "#ff00ff"

            //return reads name without "/1", "/2" end.
            static std::string cutReadName(seqan::BamAlignmentRecord read);

            //translate dna5 to string
            static std::string dna5ToString(seqan::Dna5 *seq, int len);
        public:
            GraphBuilder() = default;

            //fun that need to call for add conection between contigs;
            virtual void evaluate() = 0;

            //set graph, and change lib in it to new.
            // need to be coll after setLibName
            virtual void setGraph(ContigGraph *graph);

            //set libName
            void setLibName(std::string libName);
        private:
            DECL_LOGGER("GraphBuilder");
        };
    }
}

#endif //SCAFFOLDER_GRAPHBUILDER_H

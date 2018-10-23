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
            ContigGraph::Lib::Type  libType;

            std::string getLibColor() {
                return "#000000";
            }
            ContigGraph::Lib::Type getLibType() {
                return libType;
            }

            //return reads name without "/1", "/2" end.
            static std::string cutReadName(seqan::BamAlignmentRecord read);

            //translate dna5 to string
            static std::string dna5ToString(seqan::Dna5 *seq, int len);
        public:
            GraphBuilder() = default;

            //function which need to call for add conection between contigs;
            virtual void evaluate() = 0;

            //set graph, and change lib in it to new.
            // need to be coll after setLibName
            virtual void setGraph(ContigGraph *graph);

            //set libName
            void setLibName(std::string libName);
            void setLibType(ContigGraph::Lib::Type ltype) {
                libType = ltype;
            }
        private:
            DECL_LOGGER("GraphBuilder");
        };
    }
}

#endif //SCAFFOLDER_GRAPHBUILDER_H

#ifndef SCAFFOLDER_FASTGGRAPHBUILDER_H
#define SCAFFOLDER_FASTGGRAPHBUILDER_H

#include <seqan/seq_io.h>
#include "GraphBuilder.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
        /*
         * Generate contigs and build connection connection by FASTG file
         */
        class FASTGGraphBuilder : public GraphBuilder {
        private:
            std::string fastgFileName;

            std::string parseFirstEdge(std::string s);
            std::vector<std::string> parseFollowingEdges(std::string s);
        public:
            void setFASTGFile(std::string file_name);

            void evaluate() override;

            void parseFASTG();
        private:
            DECL_LOGGER("FASTGGraphBuilder");
        };
    }
}

#endif //SCAFFOLDER_FASTGGRAPHBUILDER_H

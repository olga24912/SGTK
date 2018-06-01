#ifndef SCAFFOLDER_FASTGGRAPHBUILDER_H
#define SCAFFOLDER_FASTGGRAPHBUILDER_H

#include <seqan/seq_io.h>
#include "GraphBuilder.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
        class FASTGGraphBuilder : public GraphBuilder {
        private:
            std::map<std::string, int> contigsId;
            std::vector<std::string> contigsName;

            std::string fastgFileName;
            std::string contigFileName;

            void initGraph();
            std::string parseFirstEdge(std::string s);
            std::vector<std::string> parseFollowingEdges(std::string s);
        public:
            void setFASTGFile(std::string file_name);
            void setContigFile(std::string file_name);

            void evaluate() override;

            void parseFASTG();
        };
    }
}

#endif //SCAFFOLDER_FASTGGRAPHBUILDER_H

#ifndef SCAFFOLDER_CONNECTIONGRAPHBUILDER_H
#define SCAFFOLDER_CONNECTIONGRAPHBUILDER_H

#include <seqan/seq_io.h>
#include "GraphBuilder.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
        class ConnectionGraphBuilder : public GraphBuilder {
        private:
            std::map<std::string, int> contigsId;
            std::vector<std::string> contigsName;

            std::string contigFileName;
            std::string connectionFileName;

            void initGraph();
        public:
            void setConnectionFile(std::string file_name);

            void setContigFile(std::string file_name);

            void evaluate() override;

            void parseConnection();
        };
    }
}


#endif //SCAFFOLDER_CONNECTIONGRAPHBUILDER_H

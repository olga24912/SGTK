#ifndef SCAFFOLDER_CONNECTIONGRAPHBUILDER_H
#define SCAFFOLDER_CONNECTIONGRAPHBUILDER_H

#include <seqan/seq_io.h>
#include "GraphBuilder.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"

namespace builder {
    namespace graph_builder {
        /*
         * Build graph by list of connection between contigs.
         */
        class ConnectionGraphBuilder : public GraphBuilder {
        private:
            std::string connectionFileName;
        public:
            void setConnectionFile(std::string file_name);

            void evaluate() override;

            void parseConnection();
        private:
            DECL_LOGGER("ConnectionGraphBuilder");
        };
    }
}


#endif //SCAFFOLDER_CONNECTIONGRAPHBUILDER_H

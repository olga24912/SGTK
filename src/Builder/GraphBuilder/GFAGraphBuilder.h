#ifndef SCAFFOLDER_GFAGRAPHBUILDER_H
#define SCAFFOLDER_GFAGRAPHBUILDER_H

#include "GraphBuilder.h"

namespace builder {
    namespace graph_builder {
        /*
         * Generate contigs, build connection and scaffold connection by GFA file
         */
        class GFAGraphBuilder : public GraphBuilder {
        private:
            std::map<std::string, int> contigsId;
            std::vector<std::string> contigsName;

            std::string gfaFileName;

            void parseGFA();
            void handleS(std::string &name, std::string &record, bool firstLib, int i);
        public:
            void setGFAFile(std::string file_name);

            void evaluate() override;
        private:
            DECL_LOGGER("GFAGraphBuilder");
        };
    }
}


#endif //SCAFFOLDER_GFAGRAPHBUILDER_H

#ifndef SCAFFOLDER_LONGGRAPHBUILDER_H
#define SCAFFOLDER_LONGGRAPHBUILDER_H

#include "GraphBuilder.h"
#include <seqan/seq_io.h>

namespace builder {
    namespace graph_builder {
        /*
         * Generate connection between contigs by long(Nanopore) reads alignments
         */
        class LongGraphBuilder : public GraphBuilder {
        protected:
            std::string fileName;

            struct AligInfo {
                std::string pacbioName;
                int len, cb, ce;
                char strand;
                std::string nodename;
                AligInfo(std::string pacbioName, int len, int cb, int ce, char strand, std::string nodename):
                        pacbioName(pacbioName), len(len), cb(cb), ce(ce), strand(strand), nodename(nodename) {
                }
                AligInfo(){}
            };

            std::vector<AligInfo> aligInfo;

            void updateEdges(std::string s);
        public:
            LongGraphBuilder() = default;

            void setFileName(const std::string &fileName);

            virtual void evaluate();

        private:
            DECL_LOGGER("LongGraphBuilder");
        };
    }
}


#endif //SCAFFOLDER_PACBIOGRAPHBUILDER_H

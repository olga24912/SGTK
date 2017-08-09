#include "RNASplit50GraphBuilder.h"

namespace builder {
    namespace graph_builder {

        ContigGraph::Lib::Type RNASplit50GraphBuilder::getLibType() {
            return ContigGraph::Lib::RNA_SPLIT_50;
        }

        std::string RNASplit50GraphBuilder::getLibColor() {
            return "#ccccff";
        }
    }
}
#include "RNASplit30GraphBuilder.h"

namespace builder {
    namespace graph_builder {

        ContigGraph::Lib::Type RNASplit30GraphBuilder::getLibType() {
            return ContigGraph::Lib::RNA_SPLIT_30;
        }

        std::string RNASplit30GraphBuilder::getLibColor() {
            return "#eebef1";
        }
    }
}
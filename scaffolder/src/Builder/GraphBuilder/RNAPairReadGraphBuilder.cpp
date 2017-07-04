#include "RNAPairReadGraphBuilder.h"

namespace builder {
    namespace graph_builder {
        std::string RNAPairReadGraphBuilder::getLibColor() {
            TRACE("getLibColor");
            int color[3] = {rand() % 100, rand() % 150, 255};
            return colorToString(color);
        }

        ContigGraph::Lib::Type RNAPairReadGraphBuilder::getLibType() {
            TRACE("getLibType");
            return ContigGraph::Lib::RNA_PAIR;
        }
    }
}

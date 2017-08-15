#include <Filter/ContigGraph/ContigGraph.h>
#include "Rule.h"

namespace filter {
    namespace scaffolder {
        bool Rule::sameCoord1(ContigGraph *graph, int e1, int e2) {
            if (graph->getEdgeCoordB1(e1) <= graph->getEdgeCoordB1(e2) &&
                graph->getEdgeCoordB1(e2) <= graph->getEdgeCoordE1(e1)) {
                return true;
            }
            if (graph->getEdgeCoordB1(e1) <= graph->getEdgeCoordE1(e2) &&
                graph->getEdgeCoordE1(e2) <= graph->getEdgeCoordE1(e1)) {
                return true;
            }
            if (graph->getEdgeCoordB1(e2) <= graph->getEdgeCoordB1(e1) &&
                graph->getEdgeCoordB1(e1) <= graph->getEdgeCoordE1(e2)) {
                return true;
            }
            if (graph->getEdgeCoordB1(e2) <= graph->getEdgeCoordE1(e1) &&
                graph->getEdgeCoordE1(e1) <= graph->getEdgeCoordE1(e2)) {
                return true;
            }
            return false;
        }

        bool Rule::sameCoord2(ContigGraph *graph, int e1, int e2) {
            if (graph->getEdgeCoordB2(e1) <= graph->getEdgeCoordB2(e2) &&
                graph->getEdgeCoordB2(e2) <= graph->getEdgeCoordE2(e1)) {
                return true;
            }
            if (graph->getEdgeCoordB2(e1) <= graph->getEdgeCoordE2(e2) &&
                graph->getEdgeCoordE2(e2) <= graph->getEdgeCoordE2(e1)) {
                return true;
            }
            if (graph->getEdgeCoordB2(e2) <= graph->getEdgeCoordB2(e1) &&
                graph->getEdgeCoordB2(e1) <= graph->getEdgeCoordE2(e2)) {
                return true;
            }
            if (graph->getEdgeCoordB2(e2) <= graph->getEdgeCoordE2(e1) &&
                graph->getEdgeCoordE2(e1) <= graph->getEdgeCoordE2(e2)) {
                return true;
            }
            return false;
        }
    }
}
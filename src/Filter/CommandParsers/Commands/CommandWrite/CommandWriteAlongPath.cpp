#include <Filter/Writers/WriteFullGraph.h>
#include <Filter/Writers/DotWriter/ChrDotWriterBuilder.h>
#include "CommandWriteAlongPath.h"

namespace filter {
    namespace commands {
        // writeAlongPath <prefixFileName> <libNum> <dist> <minRefPathSize>
        void CommandWriteAlongPath::writeGraph(std::string argv, State &state, ContigGraph &graph) {
            std::stringstream ss(argv);
            std::string fileName;
            ss >> fileName;

            writers::ChrDotWriterBuilder writerBuilder;
            writerBuilder.setFilter(&graph);
            writerBuilder.setValidator(state.validator);
            writerBuilder.setMaxVert(state.maxVert);
            writerBuilder.setMaxEdge(state.maxEdge);
            writerBuilder.setCoordFile(state.coordFile);

            writers::WriteFullGraph writer(fileName, &graph, state.validator, &writerBuilder);

            writer.write();
        }
    }
}
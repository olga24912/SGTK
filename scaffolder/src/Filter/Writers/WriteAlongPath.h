#ifndef SCAFFOLDER_WRITEALONGPATH_H
#define SCAFFOLDER_WRITEALONGPATH_H

#include "Writer.h"

namespace filter {
    namespace writers {
//write only edges and vertex in local area from some path
        class WriteAlongPath : public Writer {
        private:
            std::string fileName; //file for write the graph in dot format
            int libId; //lib id for find path using only this edges
            int dist; //area from path
            int minSize; //minLen of path for writing
        public:
            WriteAlongPath(std::string fileName, int libId, int dist, int minSize, ContigGraph *graph1,
                           FileValidator *validator,
                           DotWriterBuilder *builder);

            void write() override;

        private:
            DECL_LOGGER("WriteAlongPath");
        };
    }
}

#endif //SCAFFOLDER_WRITEALONGPATH_H

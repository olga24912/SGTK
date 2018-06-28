#ifndef SCAFFOLDER_WRITEFULLGRAPH_H
#define SCAFFOLDER_WRITEFULLGRAPH_H

#include "Writer.h"

namespace filter {
    namespace writers {
//Write full graph
        class WriteFullGraph : public Writer {
        private:
            std::string fileName; //fileName for write full graph
        public:
            WriteFullGraph(std::string fileName, ContigGraph *graph1, FileValidator *validator,
                           DotWriterBuilder *builder);

            void write() override;

        private:
            DECL_LOGGER("WriteFullGraph");
        };
    }
}

#endif //SCAFFOLDER_WRITEFULLGRAPH_H

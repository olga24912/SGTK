#ifndef SCAFFOLDER_CHRDOTWRITER_H
#define SCAFFOLDER_CHRDOTWRITER_H

#include <Filter/Writers/Searcher.h>
#include "DotWriter.h"

namespace filter {
    namespace writers {
        class ChrDotWriter : public DotWriter {
        private:
            Searcher searcher;
        public:
            void writeVertexSet(std::vector<int> vert, std::string fileName) override;

            ChrDotWriter(ContigGraph *filter, FileValidator *validator, std::string coordFile) : searcher(filter) {
                DotWriter::graph = filter;
                DotWriter::validator = validator;
                DotWriter::coordFile = coordFile;
                aligInfo.parseCoordFile(filter, coordFile);
            }

        protected:
            void writeOneEdge(int e, std::ofstream &out) override;

        public:

            void writeOneChr(const std::string chr, std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment >> verts,
                     std::string fileName);

            void writeOnePart(const std::string chrName, std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment > > chrV, std::vector<int> allVert,  std::string fileName);
        };
    }
}


#endif //SCAFFOLDER_CHRDOTWRITER_H

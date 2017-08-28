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

            void writeOneChr(const std::string chr,
                             std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment >> verts,
                             const std::string chr2,
                             std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment >> verts2,
                             std::string fileName);

            void writeOnePart(const std::string chrName,
                              std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment > > chrV,
                              std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment > > chrV2,
                              std::vector<int> allVert,  std::string fileName);

            void findVertsForCurFile(int &pos, std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment>> &verts,
                             std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment>> &verts2,
                             std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment>> &chrv,
                             std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment>> &chrv2,
                             std::vector<int> &toPrint);

            std::ofstream &writeOneStrand(std::vector<std::pair<int, statistics::InfoAboutContigsAlig::Alignment>> chrV,
                                                      std::vector<bool, std::allocator<bool>> hasOtherEdge, std::vector<int> coord,
                                                      std::vector<std::pair<int, int>> vertE, std::ofstream& out);
        };
    }
}


#endif //SCAFFOLDER_CHRDOTWRITER_H

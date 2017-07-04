#ifndef SCAFFOLDER_READSSPLITER_H
#define SCAFFOLDER_READSSPLITER_H

#include <string>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <Logger/logger.hpp>

namespace builder {
    namespace reads_splitter {
//abstract class for splitting reads
        class ReadsSplitter {
        protected:
            //split one read with name - readName, data of this read - seq,
            //name for first part will be readName/1, for second - readName/2
            //split on part len and seq.size() - len
            //out1 for first read part, out2 - for second
            void
            splitRead(std::string readName, std::string seq, int len, seqan::SeqFileOut &out1, seqan::SeqFileOut &out2);

        public:
            //take read from rnaFileName
            //first part of split read will be in resFileName1
            //second - resFileName2
            virtual void splitReads(std::string rnaFileName, std::string resFileName1, std::string resFileName2) = 0;

        private:
            DECL_LOGGER("ReadsSplitter");
        };
    }
}

#endif //SCAFFOLDER_READSSPLITER_H

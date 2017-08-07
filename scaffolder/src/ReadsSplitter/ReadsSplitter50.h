#ifndef SCAFFOLDER_READSSPLITTER50_H
#define SCAFFOLDER_READSSPLITTER50_H

#include "ReadsSplitter.h"
#include "Builder/ContigGraph/SeqanUtils.h"
#include <seqan/seq_io.h>

namespace reads_splitter {
//split unmapped reads on two equal parts
    class ReadsSplitter50 : public ReadsSplitter {
    public:
        virtual void splitReads(std::string rnaFileName,
                                std::string resFileName1, std::string resFileName2);

    private:
        DECL_LOGGER("ReadsSplitter50");
    };
}

#endif //SCAFFOLDER_READSSPLITTER50_H

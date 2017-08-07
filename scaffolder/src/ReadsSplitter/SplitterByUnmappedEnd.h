#ifndef SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H
#define SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H

#include "ReadsSplitter.h"
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

namespace reads_splitter {
//if some read mapped on ref not fully, cut unmmaped end
    class SplitterByUnmappedEnd : public ReadsSplitter {
    private:
        const int MIN_READ_LEN = 20;
    public:
        virtual void splitReads(std::string rnaFileName,
                                std::string resFileName1, std::string resFileName2);

    private:
        DECL_LOGGER("SplitterByUnmappedEnd");
    };
}

#endif //SCAFFOLDER_SPLITTERBYUNMAPPEDEND_H

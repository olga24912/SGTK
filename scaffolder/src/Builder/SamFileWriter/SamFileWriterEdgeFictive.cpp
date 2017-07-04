#include "SamFileWriterEdgeFictive.h"

namespace builder {
    namespace sam_file_writer {
        SamFileWriterEdgeFictive::SamFileWriterEdgeFictive(std::string dir) : SamFileWriteEdge() {}

        SamFileWriterEdgeFictive::SamFileWriterEdgeFictive() : SamFileWriteEdge() {}

        void SamFileWriterEdgeFictive::setFileIn(seqan::BamFileIn *in) {}

        void
        SamFileWriterEdgeFictive::writeEdge(int edgeID, seqan::BamAlignmentRecord read1,
                                            seqan::BamAlignmentRecord read2) {}
    }
}
#include "SamFileWriterEdgeFictive.h"

SamFileWriterEdgeFictive::SamFileWriterEdgeFictive(std::string dir) : SamFileWriteEdge() {}

SamFileWriterEdgeFictive::SamFileWriterEdgeFictive() : SamFileWriteEdge() {}

void SamFileWriterEdgeFictive::setFileIn(seqan::BamFileIn *in) {}

void SamFileWriterEdgeFictive::writeEdge(int edgeID, seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {}

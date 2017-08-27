#include "RuleCovering.h"

namespace filter {
    namespace scaffolder {
        using namespace commands;
        void RuleCovering::simplifyGraph(filter::contig_graph::ContigGraph *filter) {
            std::vector<int> verts = filter->getVertexList();
            for (int v : verts) {
                std::vector<int> edges = filter->getEdges(v);

                for (int e : edges) {
                    if (filter->getLibType(filter->getEdgeLib(e)) == ContigGraph::Lib::RNA_PAIR) {
                        if (isAlone(filter, e)) {
                            checkThisEdge(filter, e);
                        }
                    }
                }
            }

        }

        void RuleCovering::checkThisEdge(ContigGraph *graph, int e) {
            int libNum = graph->getEdgeLib(e);

            State::BamFiles bamFile;

            for (State::BamFiles bfiles: bamFiles) {
                if (bfiles.lib == libNum) {
                    bamFile = bfiles;
                }
            }

            double cover1 = getCover1(graph, bamFile, e);
            double cover2 = getCover2(graph, bamFile, e);

            //std::vector<InfoAboutContigsAlig::Alignment> alig1 = alig->getAlignment(graph->getEdgeFrom(e));
            //std::vector<InfoAboutContigsAlig::Alignment> alig2 = alig->getAlignment(graph->getEdgeTo(e));

            //int u = graph->getEdgeFrom(e);
            //int v = graph->getEdgeTo(e);

            //std::cerr << "e: " << e << " w=" << graph->getEdgeWeight(e) << " cover1=" << cover1 << " cover2=" << cover2 << "\n";
            //std::cerr << "u= " << u << " " << graph->getTargetName(u) << "\n";
            //for (auto al : alig1) {
            //    std::cerr << al.chrName << " " << al.coordBegin << " " << al.coordEnd << "\n";
            //}
            //std::cerr << "v= " << v << " " << graph->getTargetName(v) << "\n";
            //for (auto al: alig2) {
            //    std::cerr << al.chrName << " " << al.coordBegin << " " << al.coordEnd << "\n";
            //}
            //std::cerr << "\n\n\n\n";

            if (cover1 * 3 < cover2 || cover2 * 3 < cover1) {
                graph->delEdge(e);
            }
        }

        double RuleCovering::getCover1(ContigGraph *graph, State::BamFiles bamFile, int e) {
            int cb, ce;
            ce = graph->getEdgeCoordE1(e);
            cb = graph->getEdgeCoordB1(e);
            int v = graph->getEdgeFrom(e);

            seqan::CharString bamFileName = bamFile.bam1;
            seqan::CharString baiFileName = bamFile.bai1;
            seqan::CharString rName = graph->getTargetName(v);
            int isRev = 0;

            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                isRev = 1;
                rName = graph->getTargetName(v).substr(0, graph->getTargetName(v).size() - 4);
            }

            seqan::BamFileIn inFile;
            if (!open(inFile, toCString(bamFileName))) {
                ERROR("Could not open " << bamFile.bam1 << " for reading.\n");
                puts(seqan::toCString(bamFileName));
                return 0;
            }

            seqan::BamIndex<seqan::Bai> baiIndex;
            if (!open(baiIndex, toCString(baiFileName))) {
                ERROR("Could not read BAI index file " << bamFile.bai1 << "\n");
                return 0;
            }

            seqan::BamHeader header;
            readHeader(header, inFile);

            int rID = 0;
            if (!seqan::getIdByName(rID, seqan::contigNamesCache(seqan::context(inFile)), rName)) {
                ERROR("Reference sequence named " << graph->getTargetName(v) << " not known.\n");
                return 0;
            }

            //std::cerr << graph->getEdgeCoordB1(e) << " " << graph->getEdgeCoordE1(e) << " " << cb << " " << ce << "\n";

            int sum = getCover(cb, ce, inFile, baiIndex, rID, isRev);
            sum += getCover(graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, inFile, baiIndex, rID, 1-isRev);

            return sum * 1.0/(graph->getEdgeCoordE1(e)  - graph->getEdgeCoordB1(e));
        }

        double RuleCovering::getCover2(ContigGraph *graph, State::BamFiles bamFile, int e) {
            int cb, ce;
            cb = graph->getEdgeCoordB2(e);
            ce = graph->getEdgeCoordE2(e);
            int v = graph->getEdgeTo(e);

            seqan::CharString bamFileName = bamFile.bam2;
            seqan::CharString baiFileName = bamFile.bai2;
            seqan::CharString rName = graph->getTargetName(v);
            int isRev = 0;

            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                isRev = 1;
                rName = graph->getTargetName(v).substr(0, graph->getTargetName(v).size() - 4);
            }

            seqan::BamFileIn inFile;
            if (!open(inFile, toCString(bamFileName))) {
                ERROR("Could not open " << bamFile.bam2 << " for reading.\n");
                return 0;
            }

            seqan::BamIndex<seqan::Bai> baiIndex;
            if (!open(baiIndex, toCString(baiFileName))) {
                ERROR("Could not read BAI index file " << bamFile.bai2 << "\n");
                return 0;
            }

            seqan::BamHeader header;
            readHeader(header, inFile);

            int rID = 0;
            if (!seqan::getIdByName(rID, seqan::contigNamesCache(seqan::context(inFile)), rName)) {
                ERROR("Reference sequence named " << graph->getTargetName(v) << " not known.\n");
                return 0;
            }

            int sum = getCover(cb, ce, inFile, baiIndex, rID, 1 - isRev);
            sum += getCover(graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, inFile, baiIndex, rID, isRev);

            return sum * 1.0/(graph->getEdgeCoordE2(e)  - graph->getEdgeCoordB2(e));
        }

        int RuleCovering::getCover(int cb, int ce, seqan::BamFileIn &inFile, const seqan::BamIndex<seqan::Bai> &baiIndex,
                                   int rID, int isRev) const {
            //std::cerr << "getCover " << cb << " " << ce << "\n";
            ++ce;
            bool hasAlignments = false;
            if (!jumpToRegion(inFile, hasAlignments, rID, cb, ce, baiIndex)) {
                ERROR("Could not jump to " << cb << ":" << ce << "\n");
                return 0;
            }
            if (!hasAlignments) {
                return 0;
            }

            int sum = 0;
            seqan::BamAlignmentRecord record;
            while (!atEnd(inFile)) {
                readRecord(record, inFile);

                if (record.rID == -1 || record.rID > rID || record.beginPos >= ce) {
                    break;
                }

                if (record.beginPos < cb) {
                    continue;
                }

                if (hasFlagRC(record) != isRev) {
                    continue;
                }

              //  std::cerr << ce << " " << record.beginPos + (int) getAlignmentLengthInRef(record) << " "
              //            << cb << " " << record.beginPos << "\n";

                sum += std::min(ce, record.beginPos + (int) getAlignmentLengthInRef(record)) - std::max(cb, record.beginPos);
            }
            return sum;
        }
    }
}

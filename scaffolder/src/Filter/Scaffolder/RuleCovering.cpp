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

            std::fstream out;
            out.open("err", std::fstream::app);

            std::vector<InfoAboutContigsAlig::Alignment> alig1 = alig->getAlignment(graph->getEdgeFrom(e));
            std::vector<InfoAboutContigsAlig::Alignment> alig2 = alig->getAlignment(graph->getEdgeTo(e));

            int u = graph->getEdgeFrom(e);
            int v = graph->getEdgeTo(e);

            out << "e: " << e << " w=" << graph->getEdgeWeight(e) << " cover1=" << cover1 << " cover2=" << cover2 << "\n";
            out << "u= " << u << " " << graph->getTargetName(u) << "\n";
            int cb = graph->getEdgeCoordB1(e);
            int ce = graph->getEdgeCoordE1(e);
            out << cb << " " << ce << "; " << graph->getTargetLen(u) - ce << " " <<  graph->getTargetLen(u) - cb << "\n";
            for (auto al : alig1) {
                out << al.chrName << " " << al.coordBegin << " " << al.coordEnd << " (" << (al.coordEnd - al.coordBegin) * 1.0/graph->getTargetLen(u) << ")\n";
            }

            out << "v= " << v << " " << graph->getTargetName(v) << "\n";

            cb = graph->getEdgeCoordB2(e);
            ce = graph->getEdgeCoordE2(e);
            out << cb << " " << ce << "; " << graph->getTargetLen(v) - ce << " " <<  graph->getTargetLen(v) - cb << "\n";

            for (auto al: alig2) {
                out << al.chrName << " " << al.coordBegin << " " << al.coordEnd << " (" << (al.coordEnd - al.coordBegin) * 1.0/graph->getTargetLen(v) << ")\n";
            }
            out << "\n\n\n\n";
            out.close();

            if (cover1 * 20 < cover2 || cover2 * 20 < cover1) {
                graph->delEdge(e);
            }
        }

        double RuleCovering::getCover(ContigGraph *graph, std::string bam, std::string bai, int cb, int ce, int v, int isRev) {
            seqan::CharString bamFileName = bam;
            seqan::CharString baiFileName = bai;
            seqan::CharString rName = graph->getTargetName(v);

            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                rName = graph->getTargetName(v).substr(0, graph->getTargetName(v).size() - 4);
            }

            seqan::BamFileIn inFile;
            if (!open(inFile, toCString(bamFileName))) {
                ERROR("Could not open " << bam << " for reading.\n");
                puts(seqan::toCString(bamFileName));
                return 0;
            }

            seqan::BamIndex<seqan::Bai> baiIndex;
            if (!open(baiIndex, toCString(baiFileName))) {
                ERROR("Could not read BAI index file " << bai << "\n");
                return 0;
            }
            seqan::BamHeader header;
            readHeader(header, inFile);
            int rID = 0;
            if (!seqan::getIdByName(rID, seqan::contigNamesCache(seqan::context(inFile)), rName)) {
                ERROR("Reference sequence named " << graph->getTargetName(v) << " not known.\n");
                return 0;
            }

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

                sum += std::min(ce, record.beginPos + (int) getAlignmentLengthInRef(record)) - std::max(cb, record.beginPos);
            }

            return sum * 1.0/(ce - cb);
        }

        double RuleCovering::getCover1(ContigGraph *graph, State::BamFiles bamFile, int e) {
            int cb, ce;
            ce = graph->getEdgeCoordE1(e);
            cb = graph->getEdgeCoordB1(e);
            int v = graph->getEdgeFrom(e);
            std::string targetName = graph->getTargetName(v);

            double sum = 0;
            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                sum = getCover(graph, bamFile.bam1, bamFile.bai1, cb, ce, v, 1);
                sum += getCover(graph, bamFile.bam2, bamFile.bai2, graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, v, 1);

            } else {
                sum = getCover(graph, bamFile.bam1, bamFile.bai1, cb, ce, v, 0);
                sum += getCover(graph, bamFile.bam2, bamFile.bai2, graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, v, 0);
            }

            return sum;
        }

        double RuleCovering::getCover2(ContigGraph *graph, State::BamFiles bamFile, int e) {
            int cb, ce;
            ce = graph->getEdgeCoordE2(e);
            cb = graph->getEdgeCoordB2(e);
            int v = graph->getEdgeTo(e);
            std::string targetName = graph->getTargetName(v);

            double sum = 0;
            if (graph->getTargetName(v)[graph->getTargetName(v).size() - 1] == 'v') {
                sum = getCover(graph, bamFile.bam2, bamFile.bai2, cb, ce, v, 0);
                sum += getCover(graph, bamFile.bam1, bamFile.bai1, graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, v, 0);

            } else {
                sum = getCover(graph, bamFile.bam2, bamFile.bai2, cb, ce, v, 1);
                sum += getCover(graph, bamFile.bam1, bamFile.bai1, graph->getTargetLen(v) - ce, graph->getTargetLen(v) - cb, v, 1);
            }

            return sum;
        }
    }
}

#include "PairReadGraphBuilder.h"
#include <iostream>

namespace builder {
    namespace graph_builder {
        int PairReadGraphBuilder::pairTarget(int id) {
            TRACE("getPairTarget id=" << id);
            return id ^ 1;
        }

        void PairReadGraphBuilder::setFileName2(const std::string &fileName2) {
            DEBUG("setFileName2 fileName2" << fileName2);
            PairReadGraphBuilder::fileName2 = fileName2;
        }

        void PairReadGraphBuilder::setFileName1(const std::string &fileName1) {
            DEBUG("setFileName1 fileName1" << fileName1);
            PairReadGraphBuilder::fileName1 = fileName1;
        }

        void PairReadGraphBuilder::evaluate() {
            read1ByName.clear();
            INFO("START build graph");
            handleReads();
            INFO("finish build graph");
        }

        void PairReadGraphBuilder::readHeaderInit() {
            INFO("start initHeader");
            typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;

            seqan::BamHeader sam_header;
            readHeader(sam_header, bamFile1);
            TBamContext const &bamContext = seqan::context(bamFile1);
            size_t contig_num = seqan::length(contigNames(bamContext));

            DEBUG("Contig num=" << contig_num);

            std::string name;
            for (int i = 0; i < static_cast<int>(contig_num); ++i) {
                int length = contigLengths(bamContext)[i];
                name = std::string(seqan::toCString(contigNames(bamContext)[i]));
                graph->addVertex(2 * i, name, length);
                name += "-rev";
                graph->addVertex(2 * i + 1, name, length);
            }

            INFO("finish initHeader");
        }

        bool PairReadGraphBuilder::isUniqueMapRead(seqan::BamAlignmentRecord read) {
            TRACE("isUniqueMapRead");
            seqan::BamTagsDict tagsDict(read.tags);
            int tagId = 0;
            if (seqan::findTagKey(tagId, tagsDict, "NH")) {
                int mapCnt = 0;
                seqan::extractTagValue(mapCnt, tagsDict, tagId);
                if (mapCnt > 1) {
                    return false;
                }
            }
            return true;
        }

        std::pair<std::string, int> PairReadGraphBuilder::processOneFirstRead(seqan::BamAlignmentRecord read) {
            std::string readName = cutReadName(read);
            TRACE("processOneFirstRead readName=" << readName);

            assert(read1ByName.count(readName) == 0);
            int target = get1Target(read);
            if (target < 0 || seqan::hasFlagSecondary(read) || !isUniqueMapRead(read)) {
                return std::make_pair("", -1);
            }
            addInfoAboutRead(readName, target, read);
            return std::make_pair(readName, target);
        }

        int PairReadGraphBuilder::get1Target(const seqan::BamAlignmentRecord &read) const {
            TRACE("get1Target");

            bool isRev = hasFlagRC(read);
            int target = 2 * (read.rID);
            if (isRev) {
                target++;
            }
            return target;
        }

        int PairReadGraphBuilder::get2Target(const seqan::BamAlignmentRecord &read) const {
            TRACE("get2Target");

            bool isRev = hasFlagRC(read);
            int target = 2 * (read.rID);
            if (!isRev) {
                target++;
            }
            return target;
        }

        void PairReadGraphBuilder::addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read) {
            TRACE("addInfoAboutRead readName=" << readName << " target=" << target);
            read1ByName[readName] = read;
        }

        void PairReadGraphBuilder::addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read) {
            TRACE("addInfoAbout2Read readName=" << readName << " target=" << target);
            read2ByName[readName] = read;
        }

        std::pair<std::string, int> PairReadGraphBuilder::processOneSecondRead(seqan::BamAlignmentRecord read) {
            std::string readName = cutReadName(read);
            TRACE("processOneSecondRead readName=" << readName);

            int target = get2Target(read);

            if (target < 0 || hasFlagSecondary(read) || !isUniqueMapRead(read)) {
                return std::make_pair("", -1);
            }
            addInfoAbout2Read(readName, target, read);
            return std::make_pair(readName, target);
        }

        void PairReadGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {
            DEBUG("incEdgeWeight read1 " << std::string(seqan::toCString(read1.qName)) << " " << read1.beginPos << " "
                                         << (read1.beginPos + seqan::getAlignmentLengthInRef(read1)) << " RC=" <<  hasFlagRC(read1) <<
                                         " target " << get1Target(read1));
            DEBUG("incEdgeWeight read2 " << std::string(seqan::toCString(read1.qName)) << " " << read2.beginPos << " "
                                         << (read2.beginPos + seqan::getAlignmentLengthInRef(read2)) << " RC=" << hasFlagRC(read2) <<
                                         " target " << get2Target(read2));
            assert(seqan::isUniqueMapRead(read1));
            assert(seqan::isUniqueMapRead(read2));

            int target1 = get1Target(read1);
            std::pair<int, int> t1c = getCoord1(read1, target1);
            int target2 = get2Target(read2);
            std::pair<int, int> t2c = getCoord2(read2, target2);

            if (target1 == target2 || target1 == pairTarget(target2)) {
                return;
            }


            int e1 = changeEdges(target1, t1c, target2, t2c);
            int e2 = changeEdges(pairTarget(target2),
                                 std::make_pair(graph->getTargetLen(target2) - t2c.second, graph->getTargetLen(target2) - t2c.first),
                                 pairTarget(target1),
                                 std::make_pair(graph->getTargetLen(target1) - t1c.second, graph->getTargetLen(target1) - t1c.first));

            samFileWriter.writeEdge(e1, read1, read2);
            samFileWriter.writeEdge(e2, read2, read1);
        }

        void PairReadGraphBuilder::handleReads() {
            INFO("start handle reads");
            DEBUG("fileName1=" << fileName1 << " fileName2=" << fileName2);

            if (!open(bamFile1, fileName1.c_str())) {
                std::cerr << "could not open file";
                return;
            }
            open(bamFile2, fileName2.c_str());
            DEBUG("finish open files");

            samFileWriter.setFileIn(&bamFile1);

            seqan::BamHeader samHeader2;
            readHeader(samHeader2, bamFile2);

            if (graph->getLibNum() == 1) {
                readHeaderInit();
            } else {
                seqan::BamHeader samHeader1;
                readHeader(samHeader1, bamFile1);
            }
            DEBUG("finish read header");

            seqan::BamAlignmentRecord read1, read2;

            int cnt = 0;
            while (!atEnd(bamFile1) || !atEnd(bamFile2)) {
                TRACE("next read");
                std::pair<std::string, int> readInfo1;
                std::pair<std::string, int> readInfo2;

                if (!atEnd(bamFile1)) {
                    seqan::readRecord(read1, bamFile1);
                    TRACE("read first rec");
                    readInfo1 = processOneFirstRead(read1);
                }

                TRACE("finidh process first");
                if (!atEnd(bamFile2)) {
                    seqan::readRecord(read2, bamFile2);
                    TRACE("read sec rec");
                    readInfo2 = processOneSecondRead(read2);
                }

                TRACE("read1Name=" << readInfo1.first << " read2Name=" << readInfo2.first);

                if (readInfo2.first != "" && read1ByName.count(readInfo2.first)) {
                    incEdgeWeight(read1ByName[readInfo2.first], read2);
                    read2ByName.erase(readInfo2.first);
                }
                read1ByName.erase(readInfo2.first);

                if (readInfo1.first != "" && read2ByName.count(readInfo1.first)) {
                    incEdgeWeight(read1, read2ByName[readInfo1.first]);
                    read1ByName.erase(readInfo1.first);
                }
                read2ByName.erase(readInfo1.first);
                if (cnt % 1000 == 0) {
                    INFO("finish handle first " << cnt << "reads")
                    INFO("readInfo1 =" << readInfo1.first << " " << readInfo2.second);
                }
                ++cnt;
            }
            INFO("finish handle reads");

            close(bamFile1);
            close(bamFile2);
        }

        std::pair<int, int> PairReadGraphBuilder::getCoord1(seqan::BamAlignmentRecord read, int target) {
            if ((hasFlagRC(read))) {
                return std::make_pair(graph->getTargetLen(target) - (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)),
                                      graph->getTargetLen(target) - read.beginPos);
            } else {
                return std::make_pair(read.beginPos, (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)));
            }
        }

        std::pair<int, int> PairReadGraphBuilder::getCoord2(seqan::BamAlignmentRecord read, int target) {
            if (!(hasFlagRC(read))) {
                return std::make_pair(graph->getTargetLen(target) -
                        (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)), graph->getTargetLen(target) - read.beginPos);
            } else {
                return std::make_pair(read.beginPos, (int)(read.beginPos + seqan::getAlignmentLengthInRef(read)));
            }
        }

        int
        PairReadGraphBuilder::changeEdges(int v1, std::pair<int, int> c1, int v2, std::pair<int, int> c2) {
            std::vector<ContigGraph::Edge> edges = graph->getEdgesBetween(v1, v2);

            for (ContigGraph::Edge edge: edges) {
                if (isGoodEdgeFor1(edge, c1) && isGoodEdgeFor2(edge, c2)) {
                    std::pair<int, int> coord1 = relaxCoord(std::make_pair(edge.coordBegin1, edge.coordEnd1), c1);
                    std::pair<int, int> coord2 = relaxCoord(std::make_pair(edge.coordBegin2, edge.coordEnd2), c2);

                    graph->incEdge(edge.id, coord1, coord2);

                    return edge.id;
                }
            }

            int e = graph->addEdge(v1, v2, c1, c2);
            return e;
        }

        bool PairReadGraphBuilder::isGoodEdgeFor1(ContigGraph::Edge edge, std::pair<int, int> c) {
            if (c.first < c.second) {
                if (edge.coordBegin1 >= edge.coordEnd1) return false;

                return (edge.coordBegin1 - 200 <= c.first && c.first <= edge.coordEnd1 + 200) ||
                        (edge.coordBegin1 - 200 <= c.second && c.second <= edge.coordEnd1 + 200) ||
                        (c.first - 200 <= edge.coordBegin1 && edge.coordBegin1 <= c.second + 200) ||
                        (c.first - 200 <= edge.coordEnd1 && edge.coordEnd1 <= c.second + 200);
            } else {
                if (edge.coordBegin1 <= edge.coordEnd1) return false;

                return (edge.coordEnd1 - 200 <= c.first && c.first <= edge.coordBegin1 + 200) ||
                (edge.coordEnd1 - 200 <= c.second && c.second <= edge.coordBegin1 + 200) ||
                (c.second - 200 <= edge.coordBegin1 && edge.coordBegin1 <= c.first + 200) ||
                (c.second - 200 <= edge.coordEnd1 && edge.coordEnd1 <= c.first + 200);
            }
        }

        bool PairReadGraphBuilder::isGoodEdgeFor2(ContigGraph::Edge edge, std::pair<int, int> c) {
            if (c.first < c.second) {
                if (edge.coordBegin2 >= edge.coordEnd2) return false;

                return (edge.coordBegin2 - 200 <= c.first && c.first <= edge.coordEnd2 + 200) ||
                       (edge.coordBegin2 - 200 <= c.second && c.second <= edge.coordEnd2 + 200) ||
                       (c.first - 200 <= edge.coordBegin2 && edge.coordBegin2 <= c.second + 200) ||
                       (c.first - 200 <= edge.coordEnd2 && edge.coordEnd2 <= c.second + 200);

            } else {
                if (edge.coordBegin2 <= edge.coordEnd2) return false;

                return (edge.coordEnd2 - 200 <= c.first && c.first <= edge.coordBegin2 + 200) ||
                       (edge.coordEnd2 - 200 <= c.second && c.second <= edge.coordBegin2 + 200) ||
                       (c.second - 200 <= edge.coordBegin2 && edge.coordBegin2 <= c.first + 200) ||
                       (c.second - 200 <= edge.coordEnd2 && edge.coordEnd2 <= c.first + 200);
            }
        }

        std::pair<int, int> PairReadGraphBuilder::relaxCoord(std::pair<int, int> c1, std::pair<int, int> c2) {
            if (c2.first < c2.second) {
                return std::make_pair(std::min(c2.first, c1.first), std::max(c2.second, c1.second));
            } else {
                return std::make_pair(std::max(c2.first, c1.first), std::min(c2.second, c1.second));
            }
        }
    }
}
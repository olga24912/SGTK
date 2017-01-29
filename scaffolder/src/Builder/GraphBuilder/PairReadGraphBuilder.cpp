#include "PairReadGraphBuilder.h"
#include <iostream>

int PairReadGraphBuilder::pairTarget(int id) {
    return id ^ 1;
}

void PairReadGraphBuilder::setFileName2(const std::string &fileName2) {
    PairReadGraphBuilder::fileName2 = fileName2;
}

void PairReadGraphBuilder::setFileName1(const std::string &fileName1) {
    PairReadGraphBuilder::fileName1 = fileName1;
}

void PairReadGraphBuilder::setOneSideReadFlag(bool flag) {
    oneSideRead = flag;
}

void PairReadGraphBuilder::evaluate() {
    read1ByName.clear();
    std::cerr << "START" << std::endl;
    handleReads();
}

void PairReadGraphBuilder::readHeaderInit() {
    typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;

    seqan::BamHeader sam_header;
    readHeader(sam_header, bamFile1);
    TBamContext const &bamContext = seqan::context(bamFile1);
    size_t contig_num = seqan::length(contigNames(bamContext));

    std::string name;
    for (int i = 0; i < static_cast<int>(contig_num); ++i) {
        int length = contigLengths(bamContext)[i];
        name = std::string(seqan::toCString(contigNames(bamContext)[i]));
        graph->addVertex(2*i, name, length);
        name += "-rev";
        graph->addVertex(2*i + 1, name, length);
    }
}

bool PairReadGraphBuilder::isUniqueMapRead(seqan::BamAlignmentRecord read) {
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

std::pair<std::string,int> PairReadGraphBuilder::processOneFirstRead(seqan::BamAlignmentRecord read) {
    std::string readName = SeqanUtils::cutReadName(read);
    assert(read1ByName.count(readName) == 0);

    int target = get1Target(read);

    if (target < 0 || seqan::hasFlagSecondary(read) || !isUniqueMapRead(read)) {
        return std::make_pair("", -1);
    }
    addInfoAboutRead(readName, target, read);
    return std::make_pair(readName, target);
}

int PairReadGraphBuilder::get1Target(const seqan::BamAlignmentRecord &read) const {
    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (isRev) {
        target++;
    }
    return target;
}

int PairReadGraphBuilder::get2Target(const seqan::BamAlignmentRecord &read) const {
    bool isRev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if ((!isRev && !oneSideRead) || (isRev && oneSideRead)) {
        target++;
    }
    return target;
}

void PairReadGraphBuilder::addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read) {
    read1ByName[readName] = read;
}

void PairReadGraphBuilder::addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read) {
    read2ByName[readName] = read;
}

std::pair<std::string, int> PairReadGraphBuilder::processOneSecondRead(seqan::BamAlignmentRecord read) {
    std::string readName = SeqanUtils::cutReadName(read);

    int target = get2Target(read);

    if (target < 0 || hasFlagSecondary(read) || !isUniqueMapRead(read)) {
        return std::make_pair("", -1);
    }
    addInfoAbout2Read(readName, target, read);
    return std::make_pair(readName, target);
}

void PairReadGraphBuilder::incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2) {
    assert(seqan::isUniqueMapRead(read1));
    assert(seqan::isUniqueMapRead(read2));

    int target1 = get1Target(read1);
    int target2 = get2Target(read2);
    if (target1 == target2 || target1 == pairTarget(target2)) {
        return;
    }

    std::cerr << target1 << " " << target2 << std::endl;

    int verFID = target1, verSID = target2, verRFID = pairTarget(verFID), verRSID = pairTarget(verSID);

    int e1 = graph->incEdgeWeight(verFID, verSID);
    int e2 = graph->incEdgeWeight(verRSID, verRFID);

    samFileWriter.writeEdge(e1, read1, read2);
    samFileWriter.writeEdge(e2, read2, read1);
}

void PairReadGraphBuilder::handleReads() {
    open(bamFile1, fileName1.c_str());
    open(bamFile2, fileName2.c_str());

    samFileWriter.setFileIn(&bamFile1);

    seqan::BamHeader samHeader2;
    readHeader(samHeader2, bamFile2);

    if (graph->getLibNum() == 1) {
        readHeaderInit();
    } else {
        seqan::BamHeader samHeader1;
        readHeader(samHeader1, bamFile1);
    }

    seqan::BamAlignmentRecord read1, read2;

    while (!atEnd(bamFile1) || !atEnd(bamFile2)) {
        std::pair<std::string, int> readInfo1;
        std::pair<std::string, int> readInfo2;

        if (!atEnd(bamFile1)) {
            seqan::readRecord(read1, bamFile1);
            readInfo1 = processOneFirstRead(read1);
        }
        if (!atEnd(bamFile2)) {
            seqan::readRecord(read2, bamFile2);
            readInfo2 = processOneSecondRead(read2);
        }

        std::cerr << readInfo1.first << " " << readInfo2.first << std::endl;
        std::cerr << readInfo1.second << " " << readInfo2.second << std::endl;

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
    }

    close(bamFile1);
    close(bamFile2);
}
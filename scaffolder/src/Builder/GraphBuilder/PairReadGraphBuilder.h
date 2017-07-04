#ifndef SCAFFOLDER_PAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_PAIRREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include "Builder/SamFileWriter/SamFileWriteEdge.h"
#include "ContigGraph/SeqanUtils.h"
#include <string>
#include <unordered_map>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

//generate connection between contigs by pair reads
class PairReadGraphBuilder: public GraphBuilder {
protected:
    bool oneSideRead = false;

    std::string fileName1;
    std::string fileName2;

    seqan::BamFileIn bamFile1;
    seqan::BamFileIn bamFile2;

    std::unordered_map<std::string, seqan::BamAlignmentRecord> read1ByName;
    std::unordered_map<std::string, seqan::BamAlignmentRecord> read2ByName;

    std::pair<std::string, int> processOneFirstRead(seqan::BamAlignmentRecord read);
    std::pair<std::string, int> processOneSecondRead(seqan::BamAlignmentRecord read);
    virtual void addInfoAboutRead(std::string readName, int target, seqan::BamAlignmentRecord read);
    virtual void addInfoAbout2Read(std::string readName, int target, seqan::BamAlignmentRecord read);
    int get2Target(const seqan::BamAlignmentRecord &read) const;
    int get1Target(const seqan::BamAlignmentRecord &read) const;
    void readHeaderInit();
    bool isUniqueMapRead(seqan::BamAlignmentRecord read);

    virtual void incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);

    int pairTarget(int id);

    void handleReads();
public:
    //set sam file for first pair read alignment
    void setFileName1(const std::string &fileName1);

    //set sam file for second pair read alignment
    void setFileName2(const std::string &fileName2);

    //set info about orintation  pair read
    void setOneSideReadFlag(bool flag);

    virtual void evaluate();
private:
    DECL_LOGGER("PairReadGraphBuilder");
};

#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

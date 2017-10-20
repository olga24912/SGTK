#ifndef SCAFFOLDER_PAIRREADGRAPHBUILDER_H
#define SCAFFOLDER_PAIRREADGRAPHBUILDER_H

#include "GraphBuilder.h"
#include "Builder/SamFileWriter/SamFileWriteEdge.h"
#include "ReadsSplitter/Utils/SeqanUtils.h"
#include <string>
#include <unordered_map>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

namespace builder {
    namespace graph_builder {
//generate connection between contigs by pair reads
        class PairReadGraphBuilder : public GraphBuilder {
        protected:
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

            virtual int get2Target(const seqan::BamAlignmentRecord &read) const;

            virtual int get1Target(const seqan::BamAlignmentRecord &read) const;

            void readHeaderInit();

            bool isUniqueMapRead(seqan::BamAlignmentRecord read);

            virtual void incEdgeWeight(seqan::BamAlignmentRecord read1, seqan::BamAlignmentRecord read2);

            int pairTarget(int id);

            void handleReads();

        public:
            PairReadGraphBuilder() = default;
            //set sam file for first pair read alignment
            void setFileName1(const std::string &fileName1);

            //set sam file for second pair read alignment
            void setFileName2(const std::string &fileName2);

            virtual void evaluate();

        private:
            DECL_LOGGER("PairReadGraphBuilder");

            std::pair<int, int> getCoord1(seqan::BamAlignmentRecord read, int target1);

            std::pair<int, int> getCoord2(seqan::BamAlignmentRecord read, int target2);

            int changeEdges(int target1, std::pair<int, int> c1, int target2, std::pair<int, int> c2);

            bool isGoodEdgeFor1(ContigGraph::Edge edge, std::pair<int, int> c);

            bool isGoodEdgeFor2(ContigGraph::Edge edge, std::pair<int, int> c);

            std::pair<int, int> relaxCoord(std::pair<int, int> c1, std::pair<int, int> c2);

            int compareReads(std::pair<std::string, int> info1, std::pair<std::string, int> info2);
        };
    }
}
#endif //SCAFFOLDER_PAIRREADGRAPHBUILDER_H

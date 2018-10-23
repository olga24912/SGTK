#ifndef SCAFFOLDER_MANAGER_H
#define SCAFFOLDER_MANAGER_H

#include <string>
#include <unordered_map>
#include "Filter/CommandParsers/Commands/Command.h"
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMinContig.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMinEdgeWeight.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandSetIgnore.h>
#include <Filter/CommandParsers/Commands/CommandScaffold.h>
#include "State.h"
#include "Filter/CommandParsers/Commands/CommandUploadGraph/CommandUploadGraph.h"
#include <fstream>
#include <sstream>

namespace filter {
    namespace commands {
        using namespace contig_graph;
/*
 * Manager parse "filter_config" file, and apply command from file to graph.
 */

/*
 * Commands:
 * uploadGraph <filename>
 * uploadScaffoldsGraph <contigFile> <scaffoldFileName>
 * addInfoToGraph <infoFile> <newLibName> <colorLib>
 * addBothPath <bothPathFile> <newLibName> <colorLib>
 * minEdgeW <libNum> <weight>
 * minContig <len>
 * mergeSimplePath <contigsFileName> <outFileName>
 * setIgnore <vertexIdStart> <vertexIdFinish>
 * setIgnoreEdge <edgeId>
 * mergeLib <libNum1> <libNum2> <newLibName>
 * print <fileName>
 * exit
 * setCoordFile <coordFile>
 * setExonBlock <crdFile>
 * setBamFile <bam1FileName> <bam2FileName> <bai1FileName> <bai2FileName> <lib>
 * statisticCorrectConnection <coordFile> <libNum>
 * statisticWeight <coordFile> <libNum> <step> <mxVal>
 * statisticDif <coordFile> <libNum> <step> <mxVal>
 * statisticWeightDif <coordFile> <libNum> <stepW> <mxW> <stepD> <mxD>
 * statisticTwoLib <coordFile> <libNum1> <step1> <mxW1> <step1_2> <mxW1_2> <libNum2> <step2> <mxW2>
 * statisticTwoCompetitor <coordFile> <libNum1> <w1> <libNum2> <w2>
 * statisticCluster <coordFile>
 * statisticWrongRight <outputFile> <libForWrong> <libForRight>*
 * statisticStrand <gffFileRef> <coordFile> <gffFileContig> <outputFile>
 * histogram <lib> <step>
 * delCoord <lib>
 */
        class Manager {
        private:
            static const std::string UPLOAD_GRAPH;
            static const std::string UPLOAD_SCAFFOLDS_GRAPH;
            static const std::string ADD_INFO_TO_GRAPH;
            static const std::string ADD_BOTH_PATH;
            static const std::string MIN_EDGE_WEIGHT;
            static const std::string MIN_CONTIG_LEN;
            static const std::string MERGE_SIMPLE_PATH;
            static const std::string SET_IGNORE;
            static const std::string SET_IGNORE_EDGE;
            static const std::string MERGE_LIB;
            static const std::string PRINT;
            static const std::string EXIT;
            static const std::string SET_COORD_FILE;
            static const std::string SET_EXON_BLOCK_FILE;
            static const std::string SET_BAM_FILE;
            static const std::string STAT_CORRECT_CON;
            static const std::string STAT_WEIGHT;
            static const std::string STAT_DIF;
            static const std::string STAT_WEIGHT_DIF;
            static const std::string STAT_TWO_LIB;
            static const std::string STAT_TWO_COMP;
            static const std::string STAT_CLUST;
            static const std::string STAT_WRONG_RIGHT;
            static const std::string STAT_STRAND;
            static const std::string HISTOGRAM;
            static const std::string DEL_COORD;


            std::string configFile = "filter_config";
            State state;

            ContigGraph graph;

            std::unordered_map<std::string, Command *> commandByKeyWord;

            void readConfig();

            bool handlingRequest(std::istream &in);

        public:
            Manager();

            void main(std::string configFileName);

            ~Manager() {
                for (auto it = commandByKeyWord.begin(); it != commandByKeyWord.end(); ++it) {
                    delete (it->second);
                }
            }

        private:
            DECL_LOGGER("Manager");
        };
    }
}

#endif //SCAFFOLDER_MANAGER_H
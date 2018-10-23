#include <Filter/CommandParsers/Commands/CommandPrint.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMergeLib.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandCorrectConnectionStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandWeightStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandDifStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandWeightDifStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandTwoLibStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandHistogram.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandTwoCompStatistic.h>
#include <Filter/CommandParsers/Commands/CommandSetCoordFile.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandClusterStatistic.h>
#include <Filter/CommandParsers/Commands/CommandSetBamFile.h>
#include <Filter/CommandParsers/Commands/CommandUploadGraph/CommandUploadScaffoldsGraph.h>
#include <Filter/CommandParsers/Commands/CommandUploadGraph/CommandAddInfoToGraph.h>
#include <Filter/CommandParsers/Commands/CommandUploadGraph/CommandAddBothPath.h>
#include <Filter/CommandParsers/Commands/CommandSetExonBlockFile.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandWrongRightStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandStrandStatistic.h>
#include <Filter/CommandParsers/Commands/CommandDelCoord.h>
#include "Manager.h"

namespace filter {
    namespace commands {
        const std::string Manager::UPLOAD_GRAPH = "uploadGraph";
        const std::string Manager::UPLOAD_SCAFFOLDS_GRAPH = "uploadScaffoldsGraph";
        const std::string Manager::ADD_INFO_TO_GRAPH = "addInfoToGraph";
        const std::string Manager::ADD_BOTH_PATH = "addBothPath";
        const std::string Manager::MIN_EDGE_WEIGHT = "minEdgeW";
        const std::string Manager::MIN_CONTIG_LEN = "minContig";
        const std::string Manager::SET_IGNORE = "setIgnore";
        const std::string Manager::SET_IGNORE_EDGE = "setIgnoreEdge";
        const std::string Manager::MERGE_SIMPLE_PATH = "mergeSimplePath";
        const std::string Manager::MERGE_LIB = "mergeLib";
        const std::string Manager::PRINT = "print";
        const std::string Manager::EXIT = "exit";
        const std::string Manager::SET_COORD_FILE = "setCoordFile";
        const std::string Manager::SET_EXON_BLOCK_FILE = "setExonBlock";
        const std::string Manager::SET_BAM_FILE = "setBamFile";
        const std::string Manager::STAT_CORRECT_CON = "statisticCorrectConnection";
        const std::string Manager::STAT_WEIGHT = "statisticWeight";
        const std::string Manager::STAT_DIF = "statisticDif";
        const std::string Manager::STAT_WEIGHT_DIF = "statisticWeightDif";
        const std::string Manager::STAT_TWO_LIB = "statisticTwoLib";
        const std::string Manager::STAT_TWO_COMP = "statisticTwoCompetitor";
        const std::string Manager::STAT_CLUST = "statisticCluster";
        const std::string Manager::STAT_WRONG_RIGHT = "statisticWrongRight";
        const std::string Manager::STAT_STRAND = "statisticStrand";
        const std::string Manager::HISTOGRAM = "histogram";
        const std::string Manager::DEL_COORD = "delCoord";

        Manager::Manager() {
            graph = ContigGraph();

            commandByKeyWord[UPLOAD_GRAPH] = new CommandUploadGraph();
            commandByKeyWord[UPLOAD_SCAFFOLDS_GRAPH] = new CommandUploadScaffoldsGraph();
            commandByKeyWord[ADD_INFO_TO_GRAPH] = new CommandAddInfoToGraph();
            commandByKeyWord[ADD_BOTH_PATH] = new CommandAddBothPath();
            commandByKeyWord[MIN_CONTIG_LEN] = new CommandMinContig();
            commandByKeyWord[MIN_EDGE_WEIGHT] = new CommandMinEdgeWeight();
            commandByKeyWord[SET_IGNORE] = new CommandSetIgnore();
            commandByKeyWord[SET_IGNORE_EDGE] = new CommandSetIgnoreEdge();
            commandByKeyWord[MERGE_SIMPLE_PATH] = new CommandScaffold();
            commandByKeyWord[MERGE_LIB] = new CommandMergeLib();
            commandByKeyWord[PRINT] = new CommandPrint();
            commandByKeyWord[SET_COORD_FILE] = new CommandSetCoordFile();
            commandByKeyWord[SET_EXON_BLOCK_FILE] = new CommandSetExonBlockFile();
            commandByKeyWord[SET_BAM_FILE] = new CommandSetBamFile();
            commandByKeyWord[STAT_CORRECT_CON] = new CommandCorrectConnectionStatistic();
            commandByKeyWord[STAT_WEIGHT] = new CommandWeightStatistic();
            commandByKeyWord[STAT_DIF] = new CommandDifStatistic();
            commandByKeyWord[STAT_WEIGHT_DIF] = new CommandWeightDifStatistic();
            commandByKeyWord[STAT_TWO_LIB] = new CommandTwoLibStatistic();
            commandByKeyWord[STAT_TWO_COMP] = new CommandTwoCompStatistic();
            commandByKeyWord[STAT_CLUST] = new CommandClusterStatistic();
            commandByKeyWord[STAT_WRONG_RIGHT] = new CommandWrongRightStatistic();
            commandByKeyWord[STAT_STRAND] = new CommandStrandStatistic();
            commandByKeyWord[HISTOGRAM] = new CommandHistogram();
            commandByKeyWord[DEL_COORD] = new CommandDelCoord();
        }

        void Manager::main(std::string configFileName) {
            if (configFileName != "") {
                configFile = configFileName;
            }
            readConfig();
        }

        void Manager::readConfig() {
            std::ifstream in(configFile);

            while (handlingRequest(in));
        }

        bool Manager::handlingRequest(std::istream &in) {
            std::string keyWord;
            std::string argv;

            if (!(in >> keyWord)) {
                return false;
            }


            DEBUG("key word=" << keyWord);

            if (commandByKeyWord.count(keyWord) > 0) {
                getline(in, argv);

                DEBUG("correct keyWord=" << keyWord);
                commandByKeyWord[keyWord]->execute(argv, state, graph);
            } else if (keyWord == Manager::EXIT) {
                return false;
            }
            return true;
        }
    }
}


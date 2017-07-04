#include <Filter/CommandParsers/Commands/CommandPrint.h>
#include <Filter/Filters/FilterMergeLib.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMergeLib.h>
#include <Filter/CommandParsers/Commands/CommandFV/CommandFVWithDifInLib.h>
#include <Filter/CommandParsers/Commands/CommandSetMaxVEinOneFile.h>
#include <Filter/CommandParsers/Commands/CommandFV/CommandSetFVFork.h>
#include <Filter/CommandParsers/Commands/CommandFV/CommandSetFVOnlyFirst.h>
#include <Filter/Filters/FilterIgnoreEdge.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandCorrectConnectionStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandWeightStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandDifStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandWeightDifStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandTwoLibStatistic.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandHistogram.h>
#include <Filter/CommandParsers/Commands/CommandStatistic/CommandTwoCompStatistic.h>
#include <Filter/CommandParsers/Commands/CommandFV/CommandSetFVFewParts.h>
#include <Filter/CommandParsers/Commands/CommandDW/CommandSetBSDotWriter.h>
#include "Manager.h"

const std::string Manager::UPLOAD_GRAPH = "uploadGraph";
const std::string Manager::MIN_EDGE_WEIGHT = "minEdgeW";
const std::string Manager::MIN_CONTIG_LEN = "minContig";
const std::string Manager::SET_IGNORE = "setIgnore";
const std::string Manager::RESET_IGNORE = "resetIgnore";
const std::string Manager::SET_IGNORE_EDGE = "setIgnoreEdge";
const std::string Manager::WRITE_FULL = "writeFull";
const std::string Manager::WRITE_LOCAL = "writeLocal";
const std::string Manager::WRITE_ALL_LOCAL = "writeAllLocal";
const std::string Manager::WRITE_LOCAL_VERT_IN_SEG = "writeLocalSeg";
const std::string Manager::WRITE_BIG_COMP = "writeBig";
const std::string Manager::WRITE_LOCAL_ALONG_PATH = "writeAlongPath";
const std::string Manager::MERGE_SIMPLE_PATH = "mergeSimplePath";
const std::string Manager::MERGE_LIB = "mergeLib";
const std::string Manager::PRINT = "print";
const std::string Manager::EXIT = "exit";
const std::string Manager::SET_FV_NOT_PATH_WITH_ALL_LIB = "setFileVNotPathWithAllLib";
const std::string Manager::SET_FV_WITH_DIF_IN_LIB = "setFileVWithDifInLib";
const std::string Manager::SET_FV_FORK = "setFileVFork";
const std::string Manager::SET_FV_ONLY_FIRST = "setFileVOnlyFirst";
const std::string Manager::SET_FV_FEW_PARTS = "setFileVFewParts";
const std::string Manager::SET_BLOCK_SPLIT_DOT_WRITER = "setBlockSplitDotWriter";
const std::string Manager::SET_MAX_VE_IN_ONE_FILE = "setMaxVEinOneFile";
const std::string Manager::STAT_CORRECT_CON = "statisticCorrectConnection";
const std::string Manager::STAT_WEIGHT = "statisticWeight";
const std::string Manager::STAT_DIF = "statisticDif";
const std::string Manager::STAT_WEIGHT_DIF = "statisticWeightDif";
const std::string Manager::STAT_TWO_LIB = "statisticTwoLib";
const std::string Manager::STAT_TWO_COMP = "statisticTwoCompetitor";
const std::string Manager::HISTOGRAM = "histogram";

const std::string Manager::CONFIG_FILE = "filter_config";

Manager::Manager() {
    filter = new FilterIgnoreEdge(new FilterIgnore(new FilterMinWeight
                                                           (new FilterMergeLib(new FilterAdapter(contig_graph::ContigGraph())))));

    commandByKeyWord[UPLOAD_GRAPH] = new CommandUploadGraph();
    commandByKeyWord[MIN_CONTIG_LEN] = new CommandMinContig();
    commandByKeyWord[MIN_EDGE_WEIGHT] = new CommandMinEdgeWeight();
    commandByKeyWord[SET_IGNORE] = new CommandSetIgnore();
    commandByKeyWord[RESET_IGNORE] = new CommandResetIgnore();
    commandByKeyWord[SET_IGNORE_EDGE] = new CommandSetIgnoreEdge();
    commandByKeyWord[WRITE_FULL] = new CommandWriteFull();
    commandByKeyWord[WRITE_LOCAL] = new CommandWriteLocal();
    commandByKeyWord[WRITE_ALL_LOCAL] = new CommandWriteAllLocal();
    commandByKeyWord[WRITE_LOCAL_VERT_IN_SEG] = new CommandWriteLocalVertInSeg();
    commandByKeyWord[WRITE_BIG_COMP] = new CommandWriteBigComp();
    commandByKeyWord[WRITE_LOCAL_ALONG_PATH] = new CommandWriteAlongPath();
    commandByKeyWord[MERGE_SIMPLE_PATH] = new CommandScaffold();
    commandByKeyWord[MERGE_LIB] = new CommandMergeLib();
    commandByKeyWord[PRINT] = new CommandPrint();
    commandByKeyWord[SET_FV_NOT_PATH_WITH_ALL_LIB] = new CommandSetFVNotWithAllLib();
    commandByKeyWord[SET_FV_WITH_DIF_IN_LIB] = new CommandFVWithDifInLib();
    commandByKeyWord[SET_FV_FORK] = new CommandSetFVFork();
    commandByKeyWord[SET_FV_ONLY_FIRST] = new CommandSetFVOnlyFirst();
    commandByKeyWord[SET_FV_FEW_PARTS] = new CommandSetFVFewParts();
    commandByKeyWord[SET_BLOCK_SPLIT_DOT_WRITER] = new CommandSetBSDotWriter();
    commandByKeyWord[SET_MAX_VE_IN_ONE_FILE] = new CommandSetMaxVEinOneFile();
    commandByKeyWord[STAT_CORRECT_CON] = new CommandCorrectConnectionStatistic();
    commandByKeyWord[STAT_WEIGHT] = new CommandWeightStatistic();
    commandByKeyWord[STAT_DIF] = new CommandDifStatistic();
    commandByKeyWord[STAT_WEIGHT_DIF] = new CommandWeightDifStatistic();
    commandByKeyWord[STAT_TWO_LIB] = new CommandTwoLibStatistic();
    commandByKeyWord[STAT_TWO_COMP] = new CommandTwoCompStatistic();
    commandByKeyWord[HISTOGRAM] = new CommandHistogram();
}

void Manager::main() {
    readConfig();
    INFO("Finish read config file");
}

void Manager::readConfig() {
    std::ifstream in(CONFIG_FILE);

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
        commandByKeyWord[keyWord]->execute(argv, state, filter);
    } else if (keyWord == Manager::EXIT) {
        return false;
    } else {
        std::stringstream ss(keyWord);
        int v;
        if (ss >> v && state.name == State::LOCAL) {
            std::stringstream sout;
            sout << state.fileName << " " << v << " " << state.dist;
            CommandWriteLocal().execute(string(sout.str()), state, filter);
        }

    }

    return true;
}


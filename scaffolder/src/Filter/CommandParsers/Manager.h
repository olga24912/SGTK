#ifndef SCAFFOLDER_MANAGER_H
#define SCAFFOLDER_MANAGER_H

#include <string>
#include <unordered_map>
#include <Filter/Filters/FilterIgnore.h>
#include <Filter/Filters/FilterAdapter.h>
#include <Filter/Filters/FilterMinWeight.h>
#include <Filter/Filters/Filter.h>
#include "Filter/CommandParsers/Commands/Command.h"
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMinContig.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandMinEdgeWeight.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandSetIgnore.h>
#include <Filter/CommandParsers/Commands/CommandFilter/CommandResetIgnore.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteFull.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteLocal.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteAllLocal.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteBigComp.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteLocalVertInSeg.h>
#include <Filter/CommandParsers/Commands/CommandWrite/CommandWriteAlongPath.h>
#include <Filter/CommandParsers/Commands/CommandScaffold.h>
#include <Filter/CommandParsers/Commands/CommandFV/CommandSetFVNotPathWithAllLib.h>
#include "State.h"
#include "Commands/CommandUploadGraph.h"
#include <fstream>
#include <sstream>

/*
 * commands:
 * uploadGraph <filename>
 * minEdgeW <libNum> <weight>
 * minContig <len>
 * mergeSimplePath <contigsFileName> <outFileName>
 * setIgnore <vertexIdStart> <vertexIdFinish>
 * resetIgnore
 * setIgnoreEdge <edgeId>
 * mergeLib <libNum1> <libNum2> <newLibName>
 * print <fileName>
 * exit
 * writeFull <fileName>
 * WriteLocal <fileName> <vertexID> <dist>
 * writeAllLocal <fileName> <dist>
 * writeLocalSeg <fileName> <vertexIDStart> <vertexIDFinish> <dist>
 * writeBig <prefixFileName> <size>
 * writeAlongPath <prefixFileName> <libNum> <dist> <minRefPathSize>
 * setFileVNotPathWithAllLib
 * setFileVWithDifInLib <libNum>*
 * setFileVFork <libNum>
 * setFileVOnlyFirst <libNumPresent> <libNumNotPresent>
 * setMaxVEinOneFile <maxVertNum> <maxEdgeNum>
 * statisticCorrectConnection <coordFile> <libNum>
 * statisticWeight <coordFile> <libNum> <step> <mxVal>
 * statisticDif <coordFile> <libNum> <step> <mxVal>
 * statisticWeightDif <coordFile> <libNum> <stepW> <mxW> <stepD> <mxD>
 * statisticTwoLib <coordFile> <libNum1> <step1> <mxW1> <libNum2> <step2> <mxW2>
 * histogram <lib> <step>
 */
class Manager {
private:
    static const std::string UPLOAD_GRAPH;
    static const std::string MIN_EDGE_WEIGHT;
    static const std::string MIN_CONTIG_LEN;
    static const std::string WRITE_FULL;
    static const std::string WRITE_LOCAL;
    static const std::string WRITE_LOCAL_VERT_IN_SEG;
    static const std::string WRITE_ALL_LOCAL;
    static const std::string WRITE_BIG_COMP;
    static const std::string MERGE_SIMPLE_PATH;
    static const std::string WRITE_LOCAL_ALONG_PATH;
    static const std::string SET_IGNORE;
    static const std::string RESET_IGNORE;
    static const std::string SET_IGNORE_EDGE;
    static const std::string MERGE_LIB;
    static const std::string PRINT;
    static const std::string EXIT;
    static const std::string SET_FV_NOT_PATH_WITH_ALL_LIB;
    static const std::string SET_FV_WITH_DIF_IN_LIB;
    static const std::string SET_FV_FORK;
    static const std::string SET_FV_ONLY_FIRST;
    static const std::string SET_MAX_VE_IN_ONE_FILE;
    static const std::string STAT_CORRECT_CON;
    static const std::string STAT_WEIGHT;
    static const std::string STAT_DIF;
    static const std::string STAT_WEIGHT_DIF;
    static const std::string STAT_TWO_LIB;
    static const std::string HISTOGRAM;


    static const std::string CONFIG_FILE;
    State state;

    Filter* filter;

    std::unordered_map<std::string, Command*> commandByKeyWord;
    void readConfig();
    bool handlingRequest(std::istream &in);
public:
    Manager();

    void main();

    ~Manager() {
        delete filter;
        for (auto it = commandByKeyWord.begin(); it != commandByKeyWord.end(); ++it) {
            delete (it->second);
        }
    }
};


#endif //SCAFFOLDER_MANAGER_H
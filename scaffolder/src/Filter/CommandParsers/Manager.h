#ifndef SCAFFOLDER_MANAGER_H
#define SCAFFOLDER_MANAGER_H

#include <string>
#include <unordered_map>
#include <Filter/Filters/FilterIgnore.h>
#include <Filter/Filters/FilterAdapter.h>
#include <Filter/Filters/FilterMinWeight.h>
#include <Filter/Filters/Filter.h>
#include "Filter/CommandParsers/Commands/Command.h"
#include <Filter/CommandParsers/Commands/CommandMinContig.h>
#include <Filter/CommandParsers/Commands/CommandMinEdgeWeight.h>
#include <Filter/CommandParsers/Commands/CommandSetIgnore.h>
#include <Filter/CommandParsers/Commands/CommandResetIgnore.h>
#include <Filter/CommandParsers/Commands/CommandWriteFull.h>
#include <Filter/CommandParsers/Commands/CommandWriteLocal.h>
#include <Filter/CommandParsers/Commands/CommandWriteAllLocal.h>
#include <Filter/CommandParsers/Commands/CommandWriteBigComp.h>
#include <Filter/CommandParsers/Commands/CommandWriteLocalVertInSeg.h>
#include <Filter/CommandParsers/Commands/CommandWriteAlongPath.h>
#include <Filter/CommandParsers/Commands/CommandMergeSimplePath.h>
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
 * exit
 * writeFull <fileName>
 * WriteLocal <fileName> <vertexID> <dist>
 * writeAllLocal <fileName> <dist>
 * writeLocalSeg <fileName> <vertexIDStart> <vertexIDFinish> <dist>
 * writeBig <prefixFileName> <size>
 * writeAlongPath <prefixFileName> <libNum> <dist> <minRefPathSize>
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
    static const std::string EXIT;

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
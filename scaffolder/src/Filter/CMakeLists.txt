include_directories(.)

set(SOURCE_FILES
         Scaffolder/ScafSimplePath.cpp  Scaffolder/ScafSimplePath.h
         Tools/Utils.cpp  Tools/Utils.h
         Tools/SystemTools.cpp  Tools/SystemTools.h
         Filters/Filter.h  Filters/FilterAdapter.cpp
         Filters/FilterAdapter.h  Filters/FilterMinWeight.cpp
         Filters/FilterMinWeight.h  Filters/Filter.cpp
         Filters/FilterIgnore.cpp  Filters/FilterIgnore.h
         Filters/Query.h
         Writers/DotWriter.cpp  Writers/DotWriter.h
         Writers/Searcher.cpp  Writers/Searcher.h
         Writers/Writer.h
         Writers/WriteFullGraph.cpp  Writers/WriteFullGraph.h
         Writers/WriteLocal.cpp  Writers/WriteLocal.h
         Writers/WriteBigComponent.cpp  Writers/WriteBigComponent.h
         Writers/WriteAlongPath.cpp  Writers/WriteAlongPath.h
         CommandParsers/Manager.cpp  CommandParsers/Manager.h
         CommandParsers/Commands/Command.h
         CommandParsers/State.h  CommandParsers/Commands/CommandUploadGraph.cpp
         CommandParsers/Commands/CommandUploadGraph.h
         CommandParsers/Commands/CommandMinContig.cpp  CommandParsers/Commands/CommandMinContig.h
         CommandParsers/Commands/CommandMinEdgeWeight.cpp  CommandParsers/Commands/CommandMinEdgeWeight.h
         CommandParsers/Commands/CommandSetIgnore.cpp  CommandParsers/Commands/CommandSetIgnore.h
         CommandParsers/Commands/CommandResetIgnore.cpp  CommandParsers/Commands/CommandResetIgnore.h
         CommandParsers/Commands/CommandWriteFull.cpp  CommandParsers/Commands/CommandWriteFull.h
         CommandParsers/Commands/CommandWrite.cpp  CommandParsers/Commands/CommandWrite.h
         CommandParsers/Commands/CommandWriteLocal.cpp  CommandParsers/Commands/CommandWriteLocal.h
         CommandParsers/Commands/CommandWriteAllLocal.cpp  CommandParsers/Commands/CommandWriteAllLocal.h
         CommandParsers/Commands/CommandWriteLocalVertInSeg.cpp  CommandParsers/Commands/CommandWriteLocalVertInSeg.h
         CommandParsers/Commands/CommandWriteBigComp.cpp  CommandParsers/Commands/CommandWriteBigComp.h
         CommandParsers/Commands/CommandWriteAlongPath.cpp  CommandParsers/Commands/CommandWriteAlongPath.h
         CommandParsers/Commands/CommandMergeSimplePath.cpp  CommandParsers/Commands/CommandMergeSimplePath.h)

add_library(filter_lib ${SOURCE_FILES})
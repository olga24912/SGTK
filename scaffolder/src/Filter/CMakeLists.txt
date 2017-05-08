include_directories(.)

set(SOURCE_FILES
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
         Writers/GraphSplitter.cpp Writers/GraphSplitter.h
         CommandParsers/Manager.cpp  CommandParsers/Manager.h
         CommandParsers/Commands/Command.h
         CommandParsers/State.h  CommandParsers/Commands/CommandUploadGraph.cpp
         CommandParsers/Commands/CommandUploadGraph.h
        CommandParsers/Commands/CommandFilter/CommandMinContig.cpp CommandParsers/Commands/CommandFilter/CommandMinContig.h
        CommandParsers/Commands/CommandFilter/CommandMinEdgeWeight.cpp CommandParsers/Commands/CommandFilter/CommandMinEdgeWeight.h
        CommandParsers/Commands/CommandFilter/CommandSetIgnore.cpp CommandParsers/Commands/CommandFilter/CommandSetIgnore.h
        CommandParsers/Commands/CommandFilter/CommandResetIgnore.cpp CommandParsers/Commands/CommandFilter/CommandResetIgnore.h
        CommandParsers/Commands/CommandWrite/CommandWriteFull.cpp CommandParsers/Commands/CommandWrite/CommandWriteFull.h
        CommandParsers/Commands/CommandWrite/CommandWrite.cpp CommandParsers/Commands/CommandWrite/CommandWrite.h
        CommandParsers/Commands/CommandWrite/CommandWriteLocal.cpp CommandParsers/Commands/CommandWrite/CommandWriteLocal.h
        CommandParsers/Commands/CommandWrite/CommandWriteAllLocal.cpp CommandParsers/Commands/CommandWrite/CommandWriteAllLocal.h
        CommandParsers/Commands/CommandWrite/CommandWriteLocalVertInSeg.cpp CommandParsers/Commands/CommandWrite/CommandWriteLocalVertInSeg.h
        CommandParsers/Commands/CommandWrite/CommandWriteBigComp.cpp CommandParsers/Commands/CommandWrite/CommandWriteBigComp.h
        CommandParsers/Commands/CommandWrite/CommandWriteAlongPath.cpp CommandParsers/Commands/CommandWrite/CommandWriteAlongPath.h
        CommandParsers/Commands/CommandScaffold.cpp CommandParsers/Commands/CommandScaffold.h
        CommandParsers/Commands/CommandPrint.cpp CommandParsers/Commands/CommandPrint.h
        Filters/FilterMergeLib.cpp Filters/FilterMergeLib.h
        CommandParsers/Commands/CommandFilter/CommandMergeLib.cpp
        CommandParsers/Commands/CommandFilter/CommandMergeLib.h
        Writers/FileValidator/FileValidator.h Writers/FileValidator/ValidatorNotPathWithAllLib.cpp
        Writers/FileValidator/ValidatorNotPathWithAllLib.h
        CommandParsers/Commands/CommandFV/CommandSetFV.h
        CommandParsers/Commands/CommandFV/CommandSetFVNotPathWithAllLib.h
        Writers/FileValidator/ValidatorWithDifInLib.cpp Writers/FileValidator/ValidatorWithDifInLib.h
        CommandParsers/Commands/CommandFV/CommandFVWithDifInLib.cpp
        CommandParsers/Commands/CommandFV/CommandFVWithDifInLib.h
        CommandParsers/Commands/CommandSetMaxVEinOneFile.h
        Scaffolder/Scaffolds.cpp Scaffolder/Scaffolds.h
        Scaffolder/ScaffolderPipeline.cpp Scaffolder/ScaffolderPipeline.h Writers/FileValidator/ValidatorFork.cpp Writers/FileValidator/ValidatorFork.h CommandParsers/Commands/CommandFV/CommandSetFVFork.h CommandParsers/Commands/CommandFV/CommandSetFVOnlyFirst.h Writers/FileValidator/ValidatorOnlyFirst.cpp Writers/FileValidator/ValidatorOnlyFirst.h Filters/FilterIgnoreEdge.cpp Filters/FilterIgnoreEdge.h CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.cpp CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.h Scaffolder/ScaffoldStrategy.cpp Scaffolder/ScaffoldStrategy.h Scaffolder/ScaffoldStrategyOneLine.cpp Scaffolder/ScaffoldStrategyOneLine.h Scaffolder/ScaffoldStrategyLeaveConnection.cpp Scaffolder/ScaffoldStrategyLeaveConnection.h Scaffolder/ScaffoldStrategyUniqueConnection.cpp Scaffolder/ScaffoldStrategyUniqueConnection.h Statistics/InfoAboutContigsAlig.cpp Statistics/InfoAboutContigsAlig.h)

add_library(filter_lib ${SOURCE_FILES})
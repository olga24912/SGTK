include_directories(.)

set(SOURCE_FILES
        Tools/SystemTools.cpp  Tools/SystemTools.h
        Writers/DotWriter/DotWriter.cpp Writers/DotWriter/DotWriter.h
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
        CommandParsers/Commands/CommandFilter/CommandMergeLib.cpp
        CommandParsers/Commands/CommandFilter/CommandMergeLib.h
        CommandParsers/Commands/CommandFV/CommandSetFVNotSimplePath.h
        Writers/FileValidator/FileValidator.h Writers/FileValidator/ValidatorNotPathWithAllLib.cpp
        Writers/FileValidator/ValidatorNotPathWithAllLib.h
        CommandParsers/Commands/CommandFV/CommandSetFV.h
        CommandParsers/Commands/CommandFV/CommandSetFVNotPathWithAllLib.h
        Writers/FileValidator/ValidatorWithDifInLib.cpp Writers/FileValidator/ValidatorWithDifInLib.h
        CommandParsers/Commands/CommandFV/CommandFVWithDifInLib.cpp
        CommandParsers/Commands/CommandFV/CommandFVWithDifInLib.h
        CommandParsers/Commands/CommandSetMaxVEinOneFile.h
        Scaffolder/Scaffolds.cpp Scaffolder/Scaffolds.h
        Scaffolder/ScaffolderPipeline.cpp Scaffolder/ScaffolderPipeline.h
        Writers/FileValidator/ValidatorFork.cpp Writers/FileValidator/ValidatorFork.h
        Writers/FileValidator/ValidatorNotSimplePath.cpp Writers/FileValidator/ValidatorNotSimplePath.h
        CommandParsers/Commands/CommandFV/CommandSetFVFork.h CommandParsers/Commands/CommandFV/CommandSetFVOnlyFirst.h
        Writers/FileValidator/ValidatorOnlyFirst.cpp Writers/FileValidator/ValidatorOnlyFirst.h
        CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.cpp
        CommandParsers/Commands/CommandFilter/CommandSetIgnoreEdge.h Scaffolder/ScaffoldStrategy.cpp
        Scaffolder/ScaffoldStrategy.h
        Scaffolder/ScaffoldStrategyUniqueConnection.cpp Scaffolder/ScaffoldStrategyUniqueConnection.h
        Statistics/InfoAboutContigsAlig.cpp Statistics/InfoAboutContigsAlig.h Statistics/CntRightConnectionStatistic.cpp
        Statistics/CntRightConnectionStatistic.h CommandParsers/Commands/CommandStatistic/CommandCorrectConnectionStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandCorrectConnectionStatistic.h Statistics/Statistic.cpp Statistics/Statistic.h
        Statistics/WeightStatistic.cpp Statistics/WeightStatistic.h CommandParsers/Commands/CommandStatistic/CommandWeightStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandWeightStatistic.h Statistics/DifStatistic.cpp
        Statistics/DifStatistic.h CommandParsers/Commands/CommandStatistic/CommandDifStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandDifStatistic.h Statistics/WeightDifStatistic.cpp
        Statistics/WeightDifStatistic.h CommandParsers/Commands/CommandStatistic/CommandWeightDifStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandWeightDifStatistic.h Statistics/TwoLibStatistic.cpp
        Statistics/TwoLibStatistic.h CommandParsers/Commands/CommandStatistic/CommandTwoLibStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandTwoLibStatistic.h Statistics/WeightHistogram.cpp
        Statistics/WeightHistogram.h CommandParsers/Commands/CommandStatistic/CommandHistogram.cpp
        CommandParsers/Commands/CommandStatistic/CommandHistogram.h CommandParsers/Commands/CommandStatistic/CommandTwoCompStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandTwoCompStatistic.h Statistics/TwoCompStatistic.cpp
        Statistics/TwoCompStatistic.h Writers/FileValidator/ValidatorFewParts.cpp
        Writers/FileValidator/ValidatorFewParts.h CommandParsers/Commands/CommandFV/CommandSetFVFewParts.h
        CommandParsers/Commands/CommandDW/CommandSetBSDotWriter.h Writers/DotWriter/DotWriterBuilder.h
        Writers/DotWriter/BlockSplitDotWriterBuilder.h Writers/DotWriter/BlockSplitDotWriter.cpp
        Writers/DotWriter/BlockSplitDotWriter.h Scaffolder/Rule.cpp Scaffolder/Rule.h
        Scaffolder/RuleBigDifInWeight.cpp Scaffolder/RuleBigDifInWeight.h
        Scaffolder/RuleDelCycle.cpp Scaffolder/RuleDelCycle.h
        Scaffolder/RuleInOneLine.cpp Scaffolder/RuleInOneLine.h
        Writers/FileValidator/ValidatorError.cpp Writers/FileValidator/ValidatorError.h
        CommandParsers/Commands/CommandFV/CommandSetFVError.h Scaffolder/RuleBigDeg.cpp
        Scaffolder/RuleBigDeg.h CommandParsers/Commands/CommandSetCoordFile.h Statistics/ClusterStatistic.cpp
        Statistics/ClusterStatistic.h CommandParsers/Commands/CommandStatistic/CommandClusterStatistic.cpp
        CommandParsers/Commands/CommandStatistic/CommandClusterStatistic.h
        Scaffolder/RuleDelSmallCycle.cpp Scaffolder/RuleDelSmallCycle.h
        ContigGraph/ContigGraph.h ContigGraph/ContigGraph.cpp
        ContigGraph/SeqanUtils.h ContigGraph/SeqanUtils.cpp Scaffolder/RuleDelSmallEdges.cpp Scaffolder/RuleDelSmallEdges.h)

add_library(filter_lib ${SOURCE_FILES})
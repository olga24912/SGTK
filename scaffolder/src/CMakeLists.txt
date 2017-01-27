include_directories(.)

set(SOURCE_FILES Builder/GraphControl.cpp Builder/GraphControl.h
        ContigGraph/ContigGraph.cpp ContigGraph/ContigGraph.h
        Builder/GraphBuilder/GraphBuilder.cpp Builder/GraphBuilder/GraphBuilder.h
        Builder/GraphBuilder/RNASplitReadGraphBuilder.cpp Builder/GraphBuilder/RNASplitReadGraphBuilder.h
        Builder/GraphBuilder/PairReadGraphBuilder.cpp Builder/GraphBuilder/PairReadGraphBuilder.h
        Builder/GraphBuilder/DNAPairReadGraphBuilder.cpp Builder/GraphBuilder/DNAPairReadGraphBuilder.h
        Builder/GraphBuilder/RNAPairReadGraphBuilder.cpp Builder/GraphBuilder/RNAPairReadGraphBuilder.h
        Builder/Tools/SystemAlignmentTools.cpp Builder/Tools/SystemAlignmentTools.h
        Builder/ReadsSplitter/ReadsSplitter.cpp Builder/ReadsSplitter/ReadsSplitter.h
        Builder/Tools/SeqanUtils.cpp Builder/Tools/SeqanUtils.h Builder/ExternalLib/optionparser.h
        Builder/ReadsSplitter/ReadsSplitter50.cpp Builder/ReadsSplitter/ReadsSplitter50.h
        Builder/ReadsSplitter/SplitterByUnmappedEnd.cpp Builder/ReadsSplitter/SplitterByUnmappedEnd.h
        Builder/SamFileWriter/SamFileWriteEdge.cpp Builder/SamFileWriter/SamFileWriteEdge.h
        Builder/SamFileWriter/ContigMerger.cpp Builder/SamFileWriter/ContigMerger.h
        Builder/GraphBuilder/ReferenceGraphBuilder.cpp Builder/GraphBuilder/ReferenceGraphBuilder.h
        Filter/Scaffolder/ScafSimplePath.cpp Filter/Scaffolder/ScafSimplePath.h
        Filter/Tools/Utils.cpp Filter/Tools/Utils.h
        Filter/Tools/SystemTools.cpp Filter/Tools/SystemTools.h
        Filter/Filters/Filter.h Filter/Filters/FilterAdapter.cpp
        Filter/Filters/FilterAdapter.h Filter/Filters/FilterMinWeight.cpp
        Filter/Filters/FilterMinWeight.h Filter/Filters/Filter.cpp
        Filter/Filters/FilterIgnore.cpp Filter/Filters/FilterIgnore.h
        Filter/Writers/DotWriter.cpp Filter/Writers/DotWriter.h
        Filter/Writers/Searcher.cpp Filter/Writers/Searcher.h
        Filter/Writers/Writer.h
        Filter/Writers/WriteFullGraph.cpp Filter/Writers/WriteFullGraph.h
        Filter/Writers/WriteLocal.cpp Filter/Writers/WriteLocal.h
        Filter/CommandParsers/Manager.cpp Filter/CommandParsers/Manager.h
        Filter/CommandParsers/Commands/Command.h
        Filter/CommandParsers/State.h Filter/CommandParsers/Commands/CommandUploadGraph.cpp
        Filter/CommandParsers/Commands/CommandUploadGraph.h
        Filter/Filters/Query.h
        Filter/CommandParsers/Commands/CommandMinContig.cpp Filter/CommandParsers/Commands/CommandMinContig.h
        Filter/CommandParsers/Commands/CommandMinEdgeWeight.cpp Filter/CommandParsers/Commands/CommandMinEdgeWeight.h
        Filter/CommandParsers/Commands/CommandSetIgnore.cpp Filter/CommandParsers/Commands/CommandSetIgnore.h
        Filter/CommandParsers/Commands/CommandResetIgnore.cpp Filter/CommandParsers/Commands/CommandResetIgnore.h
        Filter/CommandParsers/Commands/CommandWriteFull.cpp Filter/CommandParsers/Commands/CommandWriteFull.h
        Filter/CommandParsers/Commands/CommandWrite.cpp Filter/CommandParsers/Commands/CommandWrite.h
        Filter/CommandParsers/Commands/CommandWriteLocal.cpp Filter/CommandParsers/Commands/CommandWriteLocal.h
        Filter/CommandParsers/Commands/CommandWriteAllLocal.cpp Filter/CommandParsers/Commands/CommandWriteAllLocal.h
        Filter/CommandParsers/Commands/CommandWriteLocalVertInSeg.cpp Filter/CommandParsers/Commands/CommandWriteLocalVertInSeg.h
        Filter/CommandParsers/Commands/CommandWriteBigComp.cpp Filter/CommandParsers/Commands/CommandWriteBigComp.h
        Filter/CommandParsers/Commands/CommandWriteAlongPath.cpp Filter/CommandParsers/Commands/CommandWriteAlongPath.h
        Filter/Writers/WriteBigComponent.cpp Filter/Writers/WriteBigComponent.h
        Filter/Writers/WriteAlongPath.cpp Filter/Writers/WriteAlongPath.h Filter/CommandParsers/Commands/CommandMergeSimplePath.cpp Filter/CommandParsers/Commands/CommandMergeSimplePath.h)

add_library(scaffolder_lib ${SOURCE_FILES})
/*#include <Builder/Tools/SystemAlignmentTools.h>
#include "Builder/ReadsSplitter/ReadsSplitter.h"
#include "Builder/ReadsSplitter/ReadsSplitter50.h"
#include "Builder/ReadsSplitter/SplitterByUnmappedEnd.h"
#include "Builder/GraphBuilder/RNAPairReadGraphBuilder.h"
#include "Builder/GraphBuilder/RNASplitReadGraphBuilder.h"
#include "GenTests/RnaPairReadsTest.h"
#include "GenTests/SplitReadsTest.h"
#include "GenTests/DnaPairReadsTest.h"
#include "gtest/gtest.h"

class ReadsSplitterTest : public  ::testing::Test {
protected:
/*    std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

    std::string tmpRefFileName = "/tmp/ref1.fasta";
    std::string tmpReadsFileName = "/tmp/reads.fasta";
    std::string tmpSamFileName = "/tmp/rna.sam";
    std::string tmpUnmappedFileName = "/tmp/Unmapped.fasta";

    ReadsSplitter *splitter;
    void SetUp() {
        SplitReadsTest gen;
        gen.genTest(sourceFileName, tmpRefFileName, tmpReadsFileName, 100);

        SystemAlignmentTools st;

        st.alignmentRNA(tmpRefFileName, tmpReadsFileName, "rna.sam", "/tmp");
        std::string command = "mv /tmp/Unmapped.out.mate1 " + tmpUnmappedFileName;
        system(command.c_str());
    }

    void checkOneShortLongRead(std::string name1, std::string name2, std::string seq1, std::string seq2, int i, int len) {
        std::stringstream ss1;
        ss1 << ">read" << i << "/1";
        ASSERT_EQ(std::string(ss1.str()), name1);

        std::stringstream ss2;
        ss2 << ">read" << i << "/2";
        ASSERT_EQ(std::string(ss2.str()), name2);

        ASSERT_EQ(seq1.size(), len);
        ASSERT_EQ(seq2.size() + seq1.size(), 100);
    }
};

TEST_F(ReadsSplitterTest, testSplitter50) {
/*    splitter = new ReadsSplitter50();
    std::string fileName1 = "/tmp/read1.fasta";
    std::string fileName2 = "/tmp/read2.fasta";

    splitter->splitReads(tmpUnmappedFileName, fileName1, fileName2);

    std::ifstream uin(tmpUnmappedFileName);
    std::ifstream f1in(fileName1);
    std::ifstream f2in(fileName2);

    std::string uname, name1, name2;
    std::string useq, seq1, seq2;

    while(uin >> uname) {
        f1in >> name1;
        f2in >> name2;
        uin >> useq;
        f1in >> seq1;
        f2in >> seq2;

        ASSERT_EQ(seq1.size(), useq.size()/2);
        ASSERT_EQ(seq2.size(), useq.size()/2);

        ASSERT_EQ(name1, uname + "/1");
        ASSERT_EQ(name2, uname + "/2");

        ASSERT_EQ(seq1 + seq2, useq);
    }

    delete splitter;
}

TEST_F(ReadsSplitterTest, testShortLong) {
    /*splitter = new SplitterByUnmappedEnd;
    std::string fileName1 = "/tmp/read1.fasta";
    std::string fileName2 = "/tmp/read2.fasta";

    splitter->splitReads(tmpSamFileName, fileName1, fileName2);

    std::ifstream f1in(fileName1);
    std::ifstream f2in(fileName2);

    std::string name1, name2;
    std::string seq1, seq2;

    for (int i = 1421, len = 79; i < 1434; ++i, --len) {
        f1in >> name1;
        f2in >> name2;
        f1in >> seq1;
        std::string x;
        if (len > 70) {
            f1in >> x;
            seq1 += x;
        }
        f2in >> seq2;

        checkOneShortLongRead(name1, name2, seq1, seq2, i, len);
    }


    for (int i = 1465, len = 69; i < 1476; ++i, ++len) {
        f1in >> name1;
        f2in >> name2;
        f1in >> seq1;
        f2in >> seq2;
        std::string x;
        if (len > 70) {
            f2in >> x;
            seq2 += x;
        }


        checkOneShortLongRead(name1, name2, seq1, seq2, i, 100 - len);
    }

    delete splitter;
}

class GraphBuilderTest : public ::testing::Test {
protected:
    std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

    std::string tmpRefFileName = "/tmp/ref.fasta";
    std::string tmpRead1FileName = "/tmp/read1.fasta";
    std::string tmpSam1FileName = "/tmp/rna1.sam";
    std::string tmpRead2FileName = "/tmp/read2.fasta";
    std::string tmpSam2FileName = "/tmp/rna2.sam";

    GraphBuilder *builder;
    ContigGraph *graph;

    void genCheckGraph() {
        SystemAlignmentTools st;
        st.alignmentRNA(tmpRefFileName, tmpRead1FileName, "rna1.sam", "/tmp");
        st.alignmentRNA(tmpRefFileName, tmpRead2FileName, "rna2.sam", "/tmp");

        graph = new ContigGraph;

        builder->setLibName("testLib", "/tmp");
        (dynamic_cast<PairReadGraphBuilder*> (builder))->setFileName1(tmpSam1FileName);
        (dynamic_cast<PairReadGraphBuilder*> (builder))->setFileName2(tmpSam2FileName);
        builder->setGraph(graph);
        builder->evaluate();

        checkGraph(20);

        delete builder;
        delete graph;
    }

    void checkGraph(int edgeW) {
        std::cerr << graph->getLibNum() << std::endl;
        graph->write("/tmp/graph.gr");

        ASSERT_EQ(graph->getLibNum(), 1);
        ASSERT_EQ(graph->getLibName(0), "testLib");
        ASSERT_EQ(graph->getVertexCount(), 20);

        for (int i = 0; i < 20; ++i) {
            std::cerr << i << " " << graph->getEdges(i).size() << std::endl;
        }

        for (int i = 0; i < 17; i = (i + 2) ^ 1) {
            ASSERT_EQ(graph->getEdges(i).size(), 1);
            ASSERT_EQ(graph->getToVertex(graph->getEdges(i)[0]), (i + 2)^1);
            ASSERT_GE(graph->getEdgeWeight(graph->getEdges(i)[0]), edgeW);
        }

        for (int i = 3; i < 20; i = (i + 2)^1) {
            ASSERT_EQ(graph->getEdgesR(i).size(), 1);
            ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(i)[0]), (i - 2)^1);
        }

        for (int i = 2; i < 17; i = (i + 2) ^ 1) {
            ASSERT_EQ(graph->getEdges(i).size(), 1);
            ASSERT_EQ(graph->getToVertex(graph->getEdges(i)[0]), (i - 2)^1);
            ASSERT_GE(graph->getEdgeWeight(graph->getEdges(i)[0]), edgeW);
        }

        for (int i = 1; i < 18; i = (i + 2) ^ 1) {
            ASSERT_EQ(graph->getEdgesR(i).size(), 1);
            ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(i)[0]), (i + 2)^1);
        }


        for (int i = 0; i < 20; i += 2) {
            std::stringstream name;
            name << "contig" << i/2;

            ASSERT_EQ(1000, graph->getTargetLength(i));
            ASSERT_EQ(std::string(name.str()), graph->getTargetName(i));
        }

        for (int i = 1; i < 20; i += 2) {
            std::stringstream name;
            name << "contig" << i/2 << "-rev";

            ASSERT_EQ(1000, graph->getTargetLength(i));
            ASSERT_EQ(std::string(name.str()), graph->getTargetName(i));
        }
    }
};

TEST_F(GraphBuilderTest, testRNAPairRead) {
    RnaPairReadsTest gen;
    gen.genTest(sourceFileName, tmpRead1FileName, tmpRead2FileName, tmpRefFileName);
    //builder = new RNAPairReadGraphBuilder;

    //genCheckGraph();
}

TEST_F(GraphBuilderTest, testDNAPairRead) {
    /*DnaPairReadsTest gen;
    gen.genTest(sourceFileName, tmpRead1FileName, tmpRead2FileName, tmpRefFileName);
    builder = new DNAPairReadGraphBuilder();
    (dynamic_cast<DNAPairReadGraphBuilder*> (builder))->setDistBetweenPairReads(200);

    genCheckGraph();
}

TEST_F(GraphBuilderTest, testRefBuilder) {
    /*GenTest gen;
    gen.genTest(sourceFileName, tmpRefFileName);

    builder = new ReferenceGraphBuilder();
    graph = new ContigGraph();
    (dynamic_cast<ReferenceGraphBuilder*> (builder))->setRefFileName(sourceFileName);
    (dynamic_cast<ReferenceGraphBuilder*> (builder))->setQueryFileName(tmpRefFileName);
    builder->setLibName("testLib", "/tmp");
    builder->setGraph(graph);
    builder->evaluate();

    std::cerr << graph->getLibNum() << std::endl;
    checkGraph(1);

    delete builder;
    delete graph;
}

class SplitReadGraphBuildTest : public ::testing::Test {
protected:
    /*std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

    std::string tmpRefFileName = "/tmp/ref.fasta";
    std::string tmpReadsFileName = "/tmp/reads.fasta";

    RNASplitReadGraphBuilder* builder;
    void SetUp() {
        SplitReadsTest gen;
        gen.genTest(sourceFileName, tmpRefFileName, tmpReadsFileName, 100);
    }
};

TEST_F(SplitReadGraphBuildTest, testSplitReadBuild) {
    /*builder = new RNASplitReadGraphBuilder();
    builder->setRefFileName(tmpRefFileName);
    builder->setRnaReadFileName(tmpReadsFileName);
    builder->setLibName("testSplit", "/tmp");
    ContigGraph *graph = new ContigGraph;
    builder->setGraph(graph);
    builder->evaluate();

    graph->write("tmp/graph.gr");
    ASSERT_EQ(graph->getLibNum(), 2);
    ASSERT_EQ(graph->getVertexCount(), 4);

    ASSERT_EQ(graph->getEdges(0).size(), 2);
    ASSERT_EQ(graph->getToVertex(graph->getEdges(0)[0]), 2);
    ASSERT_EQ(graph->getToVertex(graph->getEdges(0)[1]), 2);
    ASSERT_EQ(graph->getEdges(3).size(), 2);
    ASSERT_EQ(graph->getToVertex(graph->getEdges(3)[0]), 1);
    ASSERT_EQ(graph->getToVertex(graph->getEdges(3)[1]), 1);
    ASSERT_EQ(graph->getEdgesR(2).size(), 2);
    ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(2)[0]), 0);
    ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(2)[1]), 0);
    ASSERT_EQ(graph->getEdgesR(1).size(), 2);
    ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(1)[0]), 3);
    ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(1)[1]), 3);

    delete graph;
    delete builder;
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
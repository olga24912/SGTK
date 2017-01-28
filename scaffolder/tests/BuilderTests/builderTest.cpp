#include <Builder/Tools/SystemAlignmentTools.h>
#include "Builder/ReadsSplitter/ReadsSplitter.h"
#include "Builder/ReadsSplitter/ReadsSplitter50.h"
#include "Builder/ReadsSplitter/SplitterByUnmappedEnd.h"
#include "Builder/GraphBuilder/RNAPairReadGraphBuilder.h"
#include "GenTests/RnaPairReadsTest.h"
#include "GenTests/SplitReadsTest.h"
#include "gtest/gtest.h"

class ReadsSplitterTest : public  ::testing::Test {
protected:
    std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

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
        std::string command = "mv Unmapped.out.mate1 " + tmpUnmappedFileName;
        system(command.c_str());
    }

    void checkOneShortLongRead(std::string name1, std::string name2, std::string seq1, std::string seq2, int i, int len) {
        std::stringstream ss1;
        ss1 << ">read" << i << "/1";
        ASSERT_EQ(string(ss1.str()), name1);

        std::stringstream ss2;
        ss2 << ">read" << i << "/2";
        ASSERT_EQ(string(ss2.str()), name2);

        ASSERT_EQ(seq1.size(), len);
        ASSERT_EQ(seq2.size() + seq1.size(), 100);
    }
};

TEST_F(ReadsSplitterTest, testSplitter50) {
    splitter = new ReadsSplitter50();
    std::string fileName1 = "/tmp/read1.fasta";
    std::string fileName2 = "/tmp/read2.fasta";

    splitter->splitReads(tmpUnmappedFileName, fileName1, fileName2);

    ifstream uin(tmpUnmappedFileName);
    ifstream f1in(fileName1);
    ifstream f2in(fileName2);

    string uname, name1, name2;
    string useq, seq1, seq2;

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
    splitter = new SplitterByUnmappedEnd;
    std::string fileName1 = "/tmp/read1.fasta";
    std::string fileName2 = "/tmp/read2.fasta";

    splitter->splitReads(tmpSamFileName, fileName1, fileName2);

    ifstream f1in(fileName1);
    ifstream f2in(fileName2);

    string name1, name2;
    string seq1, seq2;

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

class RNAPairReadTest : public ::testing::Test {
protected:
    std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

    std::string tmpRefFileName = "/tmp/ref.fasta";
    std::string tmpRead1FileName = "/tmp/read1.fasta";
    std::string tmpSam1FileName = "/tmp/rna1.sam";
    std::string tmpRead2FileName = "/tmp/read2.fasta";
    std::string tmpSam2FileName = "/tmp/rna2.sam";

    RNAPairReadGraphBuilder *builder;
    ContigGraph *graph;

    void SetUp() {
        RnaPairReadsTest gen;
        gen.genTest(sourceFileName, tmpRead1FileName, tmpRead2FileName, tmpRefFileName);

        SystemAlignmentTools st;
        st.alignmentRNA(tmpRefFileName, tmpRead1FileName, "rna1.sam", "/tmp");

        st.alignmentRNA(tmpRefFileName, tmpRead2FileName, "rna2.sam", "/tmp");
    }
};

TEST_F(RNAPairReadTest, testRNAPairRead) {
    builder = new RNAPairReadGraphBuilder;
    graph = new ContigGraph;

    builder->setLibName("testLib", "/tmp");
    builder->setFileName1(tmpSam1FileName);
    builder->setFileName2(tmpSam2FileName);
    builder->setGraph(graph);
    builder->evaluate();

    ASSERT_EQ(graph->getLibNum(), 1);
    ASSERT_EQ(graph->getLibName(0), "testLib");
    ASSERT_EQ(graph->getVertexCount(), 20);

    for (int i = 0; i < 17; i = (i + 2) ^ 1) {
        ASSERT_EQ(graph->getEdges(i).size(), 1);
        ASSERT_EQ(graph->getToVertex(graph->getEdges(i)[0]), (i + 2)^1);
        ASSERT_GE(graph->getEdgeWeight(graph->getEdges(i)[0]), 20);
    }

    for (int i = 3; i < 20; i = (i + 2)^1) {
        ASSERT_EQ(graph->getEdgesR(i).size(), 1);
        ASSERT_EQ(graph->getFromVertex(graph->getEdgesR(i)[0]), (i - 2)^1);
    }

    for (int i = 2; i < 17; i = (i + 2) ^ 1) {
        ASSERT_EQ(graph->getEdges(i).size(), 1);
        ASSERT_EQ(graph->getToVertex(graph->getEdges(i)[0]), (i - 2)^1);
        ASSERT_GE(graph->getEdgeWeight(graph->getEdges(i)[0]), 20);
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

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
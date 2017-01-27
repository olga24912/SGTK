#include <Builder/Tools/SystemAlignmentTools.h>
#include "Builder/ReadsSplitter/ReadsSplitter.h"
#include "Builder/ReadsSplitter/ReadsSplitter50.h"
#include "Builder/ReadsSplitter/SplitterByUnmappedEnd.h"
#include "../lib/gtest-1.7.0/include/gtest/gtest.h"
#include "GenTests/SplitReadsTest.h"

class ReadsSplitterTest : public  ::testing::Test {
protected:
    std::string sourceFileName = "../../../resources/MG1655-K12.first10K.fasta";

    std::string tmpRefFileName = "/tmp/ref.fasta";
    std::string tmpReadsFileName = "/tmp/reads.fasta";
    std::string tmpSamFileName = "/tmp/rna.sam";
    std::string tmpUnmappedFileName = "/tmp/Unmapped.fasta";

    ReadsSplitter *splitter;
    void SetUp() {
        SystemAlignmentTools st;

        st.alignmentRNA(tmpRefFileName, tmpReadsFileName, tmpSamFileName);
        std::string command = "mv Unmapped.out.mate1 " + tmpUnmappedFileName;
        system(command.c_str());

        SplitReadsTest gen;
        
        gen.genTest(sourceFileName, tmpRefFileName, tmpReadsFileName, 100);
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

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
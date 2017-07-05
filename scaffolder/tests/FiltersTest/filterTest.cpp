#include "gtest/gtest.h"
#include <Filter/Filters/FilterAdapter.h>
#include <Filter/Filters/FilterMinWeight.h>
#include <Filter/Filters/FilterIgnore.h>
#include <Filter/Writers/WriteFullGraph.h>

class FiltersTest : public  ::testing::Test {
protected:
    std::string graphFileName;

    void SetUp() {
        graphFileName = "/tmp/somefile.gr";
    }
};

TEST_F(FiltersTest, testAdapterUploadGraph) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    const int LIB_CNT = 2;
    const int VERT_CNT = 8;

    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));

    ASSERT_EQ(adapter.getLibList().size(), LIB_CNT);
    ASSERT_EQ(adapter.getVertexCount(), VERT_CNT);
}


/*TEST_F(FiltersTest, testMinWeight) {
    ContigGraph graph;
    FilterAdapter* adapter = new FilterAdapter(graph);
    FilterMinWeight filter(adapter);

    filter.processQuery(Query(Query::UPLOAD_GRAPH, graphFileName));

    check(filter, 6, 2, std::vector<int>({0, 1, 4}),
          std::vector<int>({1, 2, 0}), std::vector<int>({0, 1, 2}));

    std::cerr <<"upload" << std::endl;
    filter.processQuery(Query(Query::MIN_CONTIG_LEN, "1"));

    check(filter, 6, 2, std::vector<int>({0, 1, 2, 3, 4, 5}),
          std::vector<int>({2, 3, 1, 1, 1, 0}), std::vector<int>({1, 1, 1, 2, 2, 1}));

    std::cerr <<"min contig" << std::endl;

    filter.processQuery(Query(Query::MIN_EDGE_WEIGHT, "0 300"));

    check(filter, 6, 2, std::vector<int>({0, 1, 2, 3, 4, 5}),
          std::vector<int>({0, 2, 1, 1, 1, 0}), std::vector<int>({1, 0, 0, 2, 1, 1}));
}

TEST_F(FiltersTest, testIgnore) {
    ContigGraph graph;
    FilterAdapter* adapter = new FilterAdapter(graph);
    FilterIgnore filter(adapter);

    filter.processQuery(Query(Query::UPLOAD_GRAPH, graphFileName));

    check(filter, 6, 2, std::vector<int>({0, 1, 2, 3, 4, 5}),
          std::vector<int>({2, 3, 1, 1, 1, 0}), std::vector<int>({1, 1, 1, 2, 2, 1}));

    filter.processQuery(Query(Query::SET_IGNORE, "2 4"));

    check(filter, 6, 2, std::vector<int>({0, 1, 4, 5}),
          std::vector<int>({1, 2, 1, 0}), std::vector<int>({0, 1, 2, 1}));

    filter.processQuery(Query(Query::RESET_IGNORE, ""));


    check(filter, 6, 2, std::vector<int>({0, 1, 2, 3, 4, 5}),
          std::vector<int>({2, 3, 1, 1, 1, 0}), std::vector<int>({1, 1, 1, 2, 2, 1}));
}*/

/*
class WritersTest : public  ::testing::Test {
protected:
    Filter* filter;
    std::string graphFileName;
    std::string outFileName = "/tmp/out";

    virtual void TearDown() {
        delete filter;
    }

    void checkFiles(std::string res) {
        std::ifstream in(outFileName + "0");
        std::string out;
        std::string x;
        while (std::getline(in, x)) {
            int cntsp = 0;
            while (cntsp < (int)x.size() && x[cntsp] == ' ') ++cntsp;
            out += x.substr(cntsp) + '\n';
        }

        ASSERT_EQ(res, out);
    }

    void SetUp() {
        graphFileName = "/tmp/somefile.gr";
        std::ofstream out(graphFileName);

        out << "2\nl 0 #ff0000 lib1 RNA_PAIR\nl 1 #00ff00 lib2 RNA_SPLIT_50\n";
        out << "6\nv 0 contig0 501\nv 1 contig1 1000\nv 2 contig2 10\n";
        out << "v 3 contig3 200\nv 4 contig4 30000\nv 5 contig5 5\n";
        out << "8\n";
        out << "e 0 0 1 0 10 coord: 400 500 10 110 \n";
        out << "e 1 0 2 0 20 coord: 350 450 1 9\n";
        out << "e 2 1 4 0 4 coord: 900 950 4 54\n";
        out << "e 3 1 4 1 5000 coord: 900 950 4 54\n";
        out << "e 4 1 3 1 13 coord: 800 910 1 152\n";
        out << "e 5 2 3 0 2500 coord: 1 5 4 54\n";
        out << "e 6 3 0 1 394 coord: 150 200 40 240\n";
        out << "e 7 4 5 1 180 coord: 29800 29999 1 4\n";

        out.close();

        filter = new FilterAdapter(ContigGraph::read(graphFileName));
    }
};

TEST_F(WritersTest, testWirteFullGraph) {
    WriteFullGraph writer(outFileName, filter, new FileValidator(), new DotWriterBuilder());
    writer.write();

    std::string res = "digraph {\n"
            "\"contig0\"[label=\" contig0 id = 0\n"
            "len = 501\"];\n"
            "\"contig1\"[label=\" contig1 id = 1\n"
            "len = 1000\"];\n"
            "\"contig2\"[label=\" contig2 id = 2\n"
            "len = 10\"];\n"
            "\"contig3\"[label=\" contig3 id = 3\n"
            "len = 200\"];\n"
            "\"contig4\"[label=\" contig4 id = 4\n"
            "len = 30000\"];\n"
            "\"contig5\"[label=\" contig5 id = 5\n"
            "len = 5\"];\n"
            "\"contig1\" -> \"contig4\" [ color = \"#00ff00\", penwidth = 4, label = \"lib2\n"
            "weight = 5000\n"
            "id = 3\" ]\n"
            "\"contig2\" -> \"contig3\" [ color = \"#ff0000\", penwidth = 4, label = \"lib1\n"
            "weight = 2500\n"
            "id = 5\" ]\n"
            "\"contig3\" -> \"contig0\" [ color = \"#00ff00\", penwidth = 3, label = \"lib2\n"
            "weight = 394\n"
            "id = 6\" ]\n"
            "\"contig4\" -> \"contig5\" [ color = \"#00ff00\", penwidth = 3, label = \"lib2\n"
            "weight = 180\n"
            "id = 7\" ]\n"
            "\"contig0\" -> \"contig2\" [ color = \"#ff0000\", penwidth = 2, label = \"lib1\n"
            "weight = 20\n"
            "id = 1\" ]\n"
            "\"contig1\" -> \"contig3\" [ color = \"#00ff00\", penwidth = 2, label = \"lib2\n"
            "weight = 13\n"
            "id = 4\" ]\n"
            "\"contig0\" -> \"contig1\" [ color = \"#ff0000\", penwidth = 2, label = \"lib1\n"
            "weight = 10\n"
            "id = 0\" ]\n"
            "\"contig1\" -> \"contig4\" [ color = \"#ff0000\", penwidth = 1, label = \"lib1\n"
            "weight = 4\n"
            "id = 2\" ]\n"
            "}\n";

    checkFiles(res);
}*/


int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
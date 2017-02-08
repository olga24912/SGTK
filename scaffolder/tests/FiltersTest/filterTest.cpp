#include "gtest/gtest.h"
#include <string>
#include <iostream>
#include <Filter/Filters/FilterAdapter.h>
#include <Filter/Filters/FilterMinWeight.h>
#include <Filter/Filters/FilterIgnore.h>

class FiltersTest : public  ::testing::Test {
protected:
    ContigGraph graph;
    std::string graphFileName;

    void SetUp() {
        graphFileName = "/tmp/somefile.gr";
        std::ofstream out(graphFileName);

        out << "2\nl 0 #ff0000 lib1\nl 1 #00ff00 lib2\n";
        out << "6\nv 0 contig0 501\nv 1 contig1 1000\nv 2 contig2 10\n";
        out << "v 3 contig3 200\nv 4 contig4 30000\nv 5 contig5 5\n";
        out << "8\n";
        out << "e 0 0 1 0 10\n";
        out << "e 1 0 2 0 20\n";
        out << "e 2 1 4 0 4\n";
        out << "e 3 1 4 1 5000\n";
        out << "e 4 1 3 1 13\n";
        out << "e 5 2 3 0 2500\n";
        out << "e 6 3 0 1 394\n";
        out << "e 7 4 5 1 180\n";

        out.close();
    }

    void check(Filter& filter, int vertCount, int libCount, std::vector<int> vert,
               std::vector<int> edgesSize, std::vector<int> edgesRSize) {
        ASSERT_EQ(filter.getVertexCount(), vertCount);
        ASSERT_EQ(filter.getLibCount(), libCount);
        std::vector<int> vertn = filter.getVertexList();
        ASSERT_EQ(vertn.size(), vert.size());
        for (int i = 0; i < (int)vert.size(); ++i) {
            ASSERT_EQ(vert[i], vertn[i]);
        }

        for (int i = 0; i < vert.size(); ++i) {
            std::vector<int> edges = filter.getEdges(vertn[i]);
            ASSERT_EQ(edgesSize[i], edges.size());

            edges = filter.getEdgesR(vertn[i]);
            ASSERT_EQ(edgesRSize[i], edges.size());
        }
        std::cerr << "ok" << std::endl;
    }
};

TEST_F(FiltersTest, testAdapter) {
    ContigGraph graph;
    FilterAdapter adapter(graph);

    adapter.processQuery(Query(Query::UPLOAD_GRAPH, graphFileName));

    check(adapter, 6, 2, std::vector<int>({0, 1, 2, 3, 4, 5}),
          std::vector<int>({2, 3, 1, 1, 1, 0}), std::vector<int>({1, 1, 1, 2, 2, 1}));
}

TEST_F(FiltersTest, testMinWeight) {
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
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
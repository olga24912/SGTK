#include "gtest/gtest.h"
#include <Filter/Filters/FilterAdapter.h>
#include <Filter/Filters/FilterMinWeight.h>
#include <Filter/Filters/FilterIgnore.h>
#include <Filter/Writers/WriteFullGraph.h>
#include <Logger/log_writers.hpp>

class FiltersTest : public  ::testing::Test {
protected:
    std::string graphFileName;

    void SetUp() {
        graphFileName = "/tmp/somefile.gr";
        logging::create_console_logger("../../../log.properties");
    }

    std::string createInfoByCoord(const int* coord) {
        std::stringstream res;
        res << "\ncoord: " << coord[0] << "-" << coord[1] << "\n" << coord[2] << "-" << coord[3];
        return res.str();
    }

    std::vector<std::vector<int > > gen_gr(const int * from, int cnt_e, int cnt_v) {
        std::vector<std::vector<int> > g(cnt_v);

        for (int i = 0; i < cnt_e; ++i) {
            g[from[i]].push_back(i);
        }

        return g;
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

TEST_F(FiltersTest, testAdapterCoord) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    const int VERT_CNT = 8;
    const int COORD[10][4] = {{900, 1000, 0,   100},
                              {850, 950,  10,  110},
                              {399, 500,  11,  111},
                              {200, 300,  100, 200},
                              {301, 401,  30,  120},
                              {400, 500,  0,   100},
                              {390, 490,  50,  150},
                              {388, 488,  0,   101},
                              {301, 401,  200, 300},
                              {380, 470,  100, 200}};

    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));

    std::vector<int> vert = adapter.getVertexList();
    ASSERT_EQ(vert.size(), VERT_CNT);

    for (int v : vert) {
        std::vector<int> edgeList = adapter.getEdges(v);
        for (int e : edgeList) {
            ASSERT_EQ(createInfoByCoord(COORD[e]), adapter.getInfo(e));
            ASSERT_EQ(COORD[e][0], adapter.getEdgeCoordB1(e));
            ASSERT_EQ(COORD[e][1], adapter.getEdgeCoordE1(e));
            ASSERT_EQ(COORD[e][2], adapter.getEdgeCoordB2(e));
            ASSERT_EQ(COORD[e][3], adapter.getEdgeCoordE2(e));
        }
    }
}

TEST_F(FiltersTest, testAdapterVertex) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    const int LIB_CNT = 2;
    const int VERT_CNT = 8;

    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));

    std::vector<int> vert = adapter.getVertexList();
    ASSERT_EQ(VERT_CNT, vert.size());
    for (int i = 0; i < VERT_CNT; ++i) {
        ASSERT_EQ(i, vert[i]);
    }

    ASSERT_EQ(VERT_CNT, adapter.getVertexCount());

    for (int i = 0; i < VERT_CNT; ++i) {
        std::stringstream curName;
        curName << "node" << i;
        ASSERT_EQ(curName.str(), adapter.getTargetName(i));
    }
}

TEST_F(FiltersTest, testAdapterEdges) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    const int LIB_CNT = 2;
    const int VERT_CNT = 8;
    const int EDGE_CNT = 10;
    const int VERT_TO[EDGE_CNT] = {2, 2, 6, 4, 2, 1, 1, 3, 3, 5};
    const int VERT_FROM[EDGE_CNT] = {0, 0, 2, 2, 4, 3, 3, 7, 5, 3};
    const int EDGES_WEIGHT[EDGE_CNT] = {4, 1, 1, 2, 1, 4, 1, 1, 2, 1};
    const int EDGES_LIB[EDGE_CNT] = {0, 1, 1, 0, 1, 0, 1, 1, 0, 1};

    std::vector<std::vector<int> > gr = gen_gr(VERT_FROM, EDGE_CNT, VERT_CNT);
    std::vector<std::vector<int> > grR =  gen_gr(VERT_TO, EDGE_CNT, VERT_CNT);

    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));

    std::vector<int> vert = adapter.getVertexList();
    for (int v : vert) {
        TRACE("v = " << v);
        std::vector<int> edges = adapter.getEdges(v);
        ASSERT_EQ(gr[v], edges);
        std::vector<int> edgesR = adapter.getEdgesR(v);
        ASSERT_EQ(grR[v], edgesR);
    }

    for (int e = 0; e < EDGE_CNT; ++e) {
        ASSERT_EQ(VERT_TO[e], adapter.getEdgeTo(e));
        ASSERT_EQ(VERT_FROM[e], adapter.getEdgeFrom(e));
        ASSERT_EQ(EDGES_WEIGHT[e], adapter.getEdgeWeight(e));
        ASSERT_EQ(EDGES_LIB[e], adapter.getEdgeLib(e));
    }
}

TEST_F(FiltersTest, testAdapterLibs) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    const int LIB_CNT = 2;
    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));

    ASSERT_EQ(LIB_CNT, adapter.getLibList().size());
    std::vector<int> libs = adapter.getLibList();
    for (int i = 0; i < LIB_CNT; ++i) {
        ASSERT_EQ(i, libs[i]);
        ASSERT_EQ(adapter.getLibType(i), contig_graph::ContigGraph::Lib::RNA_PAIR);
    }

    ASSERT_EQ("#ff0000", adapter.getLibColor(0));
    ASSERT_EQ("#00ff00", adapter.getLibColor(1));
    ASSERT_EQ("lib1", adapter.getLibName(0));
    ASSERT_EQ("lib2", adapter.getLibName(1));
}

TEST_F(FiltersTest, testAdapterWrite) {
    graphFileName = "../../../resources/filterTest/simple_adapter_graph.gr";
    std::string outFile = "/tmp/out.gr";
    contig_graph::ContigGraph graph;
    filter::FilterAdapter adapter(graph);

    adapter.processQuery(filter::Query(filter::Query::UPLOAD_GRAPH, graphFileName));
    adapter.write(outFile);

    std::ifstream sin(graphFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
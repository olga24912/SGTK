#include "gtest/gtest.h"
#include <Logger/log_writers.hpp>
#include <Filter/ContigGraph/ContigGraph.h>

using namespace filter::contig_graph;

class ContigGraphTest : public  ::testing::Test {
protected:
    std::string graphFileName;

    void SetUp() {
        graphFileName = "/tmp/somefile.gr";
        logging::create_console_logger("../../../log.properties");
    }

    std::string createInfoByCoord(const int* coord) {
        std::stringstream res;
        res << "coord: " << coord[0] << "-" << coord[1] << "\n" << coord[2] << "-" << coord[3];
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

TEST_F(ContigGraphTest, testReadGraph) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
    const int LIB_CNT = 2;
    const int VERT_CNT = 8;

    ContigGraph graph = ContigGraph::read(graphFileName);

    ASSERT_EQ(graph.getLibNum(), LIB_CNT);
    ASSERT_EQ(graph.getVertexCount(), VERT_CNT);
}

TEST_F(ContigGraphTest, testEdges) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
    const int LIB_CNT = 2;
    const int VERT_CNT = 8;
    const int EDGE_CNT = 10;
    const int VERT_TO[EDGE_CNT] = {2, 2, 6, 4, 2, 1, 1, 3, 3, 5};
    const int VERT_FROM[EDGE_CNT] = {0, 0, 2, 2, 4, 3, 3, 7, 5, 3};
    const int EDGES_WEIGHT[EDGE_CNT] = {4, 1, 1, 2, 1, 4, 1, 1, 2, 1};
    const int EDGES_LIB[EDGE_CNT] = {0, 1, 1, 0, 1, 0, 1, 1, 0, 1};

    std::vector<std::vector<int> > gr = gen_gr(VERT_FROM, EDGE_CNT, VERT_CNT);
    std::vector<std::vector<int> > grR =  gen_gr(VERT_TO, EDGE_CNT, VERT_CNT);

    ContigGraph graph = ContigGraph::read(graphFileName);

    std::vector<int> vert = graph.getVertexList();
    for (int v : vert) {
        TRACE("v = " << v);
        std::vector<int> edges = graph.getEdges(v);
        ASSERT_EQ(gr[v], edges);
        std::vector<int> edgesR = graph.getEdgesR(v);
        ASSERT_EQ(grR[v], edgesR);
    }

    for (int e = 0; e < EDGE_CNT; ++e) {
        ASSERT_EQ(VERT_TO[e], graph.getEdgeTo(e));
        ASSERT_EQ(VERT_FROM[e], graph.getEdgeFrom(e));
        ASSERT_EQ(EDGES_WEIGHT[e], graph.getEdgeWeight(e));
        ASSERT_EQ(EDGES_LIB[e], graph.getEdgeLib(e));
    }
}

TEST_F(ContigGraphTest, testCoord) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
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

    ContigGraph graph = ContigGraph::read(graphFileName);

    std::vector<int> vert = graph.getVertexList();
    ASSERT_EQ(vert.size(), VERT_CNT);

    for (int v : vert) {
        std::vector<int> edgeList = graph.getEdges(v);
        for (int e : edgeList) {
            ASSERT_EQ(createInfoByCoord(COORD[e]), graph.getInfo(e));
            ASSERT_EQ(COORD[e][0], graph.getEdgeCoordB1(e));
            ASSERT_EQ(COORD[e][1], graph.getEdgeCoordE1(e));
            ASSERT_EQ(COORD[e][2], graph.getEdgeCoordB2(e));
            ASSERT_EQ(COORD[e][3], graph.getEdgeCoordE2(e));
        }
    }
}

TEST_F(ContigGraphTest, testWrite) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
    std::string outFile = "/tmp/out.gr";
    ContigGraph graph = ContigGraph::read(graphFileName);

    graph.write(outFile);

    std::ifstream sin(graphFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

TEST_F(ContigGraphTest, testVertex) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
    const int VERT_CNT = 8;
    const int LENS[VERT_CNT] = {1000, 1000, 500, 500, 501, 501, 499, 499};

    ContigGraph graph = ContigGraph::read(graphFileName);

    std::vector<int> vert = graph.getVertexList();
    ASSERT_EQ(VERT_CNT, vert.size());
    for (int i = 0; i < VERT_CNT; ++i) {
        ASSERT_EQ(i, vert[i]);
    }

    ASSERT_EQ(VERT_CNT, graph.getVertexCount());

    for (int i = 0; i < VERT_CNT; ++i) {
        std::stringstream curName;
        curName << "node" << i;
        ASSERT_EQ(curName.str(), graph.getTargetName(i));
        ASSERT_EQ(LENS[i], graph.getTargetLen(i));
        ASSERT_EQ(i, graph.getTargetId(curName.str()));
    }
}

TEST_F(ContigGraphTest, testLibs) {
    graphFileName = "../../../resources/filterTest/simple_graph.gr";
    const int LIB_CNT = 2;
    ContigGraph graph = ContigGraph::read(graphFileName);

    ASSERT_EQ(LIB_CNT, graph.getLibNum());

    std::vector<int> libs = graph.getLibList();
    for (int i = 0; i < LIB_CNT; ++i) {
        ASSERT_EQ(i, libs[i]);
        ASSERT_EQ(graph.getLibType(i), ContigGraph::Lib::RNA_PAIR);
    }

    ASSERT_EQ("#ff0000", graph.getLibColor(0));
    ASSERT_EQ("#00ff00", graph.getLibColor(1));
    ASSERT_EQ("lib1", graph.getLibName(0));
    ASSERT_EQ("lib2", graph.getLibName(1));
}

TEST_F(ContigGraphTest, testMergeGraph) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/merge1_0_lib_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);

    graph.mergeLib(0, 1, "lib01", 1, 1);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}


TEST_F(ContigGraphTest, testMergeGraphDifLib) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/merge1_2_lib_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);

    graph.mergeLib(1, 2, "lib12", 1, 1.75);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

TEST_F(ContigGraphTest, testDelEdge) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/del_edge_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);
    graph.delEdge(8);
    graph.delEdge(20);
    graph.delEdge(15);
    graph.delEdge(3);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

TEST_F(ContigGraphTest, testAddEdge) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/add_edge_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);
    graph.addEdge(9, 1, 0, 5, 0, 0, 0, 0);
    graph.addEdge(0, 8, 0, 5, 0, 0, 0, 0);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

TEST_F(ContigGraphTest, testSetWeight) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/set_weight_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);
    graph.setWeight(6, 5);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

    std::string ss, os;
    while (getline(sin, ss)) {
        getline(oin, os);
        ASSERT_EQ(ss, os);
    }

    sin.close();
    oin.close();
}

TEST_F(ContigGraphTest, testDelVertex) {
    graphFileName = "../../../resources/filterTest/merge_lib_graph.gr";
    std::string resFileName = "../../../resources/filterTest/del_vertex_graph.gr";

    ContigGraph graph = ContigGraph::read(graphFileName);
    graph.delVertex(5);

    std::string outFile = "/tmp/out.gr";
    graph.write(outFile);

    std::ifstream sin(resFileName), oin(outFile);

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
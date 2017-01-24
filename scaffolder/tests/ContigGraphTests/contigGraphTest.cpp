#include "ContigGraph/ContigGraph.h"
#include "gtest/gtest.h"
#include <string>

class ContigGraphTest : public  ::testing::Test {
protected:
    ContigGraph graph;
    std::string tmpFileName;

    void SetUp() {
        tmpFileName = "/tmp/somefile";
    }

};

TEST_F(ContigGraphTest, addNewLib) {
    ASSERT_EQ(graph.getLibNum(), 0);
    graph.newLib("lib1", "#ff0000");
    ASSERT_EQ(graph.getLibNum(), 1);
}

TEST_F(ContigGraphTest, addVertexAndEdges) {
    graph.newLib("lib1", "#ff0000");

    graph.addVertex(0, "vert1", 10);
    graph.addVertex(1, "vert1-rev", 11);
    graph.addVertex(2, "vert2", 100);
    graph.addVertex(3, "vert2-rev", 120);

    ASSERT_EQ(graph.getVertexCount(), 4);

    graph.incEdgeWeight(0, 2);

    ASSERT_EQ(graph.getTargetLength(0), 10);
    ASSERT_EQ(graph.getTargetLength(3), 120);

    graph.newLib("lib2", "#00ff00");

    graph.incEdgeWeight(0, 2);
    graph.incEdgeWeight(0, 2);
}

TEST_F(ContigGraphTest, serialization) {
    graph.newLib("lib1", "#ff0000");
    graph.addVertex(0, "vert1", 10);
    graph.addVertex(1, "vert1-rev", 11);
    graph.addVertex(2, "vert2", 100);
    graph.addVertex(3, "vert2-rev", 120);
    graph.incEdgeWeight(0, 2);
    graph.newLib("lib2", "#00ff00");
    graph.incEdgeWeight(0, 2);
    graph.incEdgeWeight(0, 2);

    graph.write(tmpFileName);

    ContigGraph graphCopy = ContigGraph::read(tmpFileName);
    ASSERT_EQ(graphCopy.getEdges(0), graph.getEdges(0));
    ASSERT_EQ(graphCopy.getEdgesR(2), graph.getEdgesR(2));
    ASSERT_EQ(graphCopy.getTargetLength(1), graph.getTargetLength(1));
    ASSERT_EQ(graphCopy.getLibNum(), graph.getLibNum());
    ASSERT_EQ(graphCopy.getToVertex(0), graph.getToVertex(0));
    ASSERT_EQ(graphCopy.getFromVertex(1), graph.getFromVertex(1));
    ASSERT_EQ(graphCopy.getVertexCount(), graph.getVertexCount());
}

TEST_F(ContigGraphTest, write) {
    graph.newLib("lib1", "#ff0000");
    graph.addVertex(0, "vert1", 10);
    graph.addVertex(1, "vert1-rev", 11);
    graph.addVertex(2, "vert2", 100);
    graph.addVertex(3, "vert2-rev", 120);
    graph.incEdgeWeight(0, 2);
    graph.newLib("lib2", "#00ff00");
    graph.incEdgeWeight(0, 2);
    graph.incEdgeWeight(0, 2);

    graph.write(tmpFileName);
    std::ifstream in(tmpFileName);
    std::string correctFile[11] = {"2", "l 0 #ff0000 lib1", "l 1 #00ff00 lib2",
                                   "4", "v 0 vert1 10", "v 1 vert1-rev 11", "v 2 vert2 100", "v 3 vert2-rev 120",
                                   "2", "e 0 0 2 0 1", "e 1 0 2 1 2"};
    for (int i = 0; i < 11; ++i) {
        std::string s;
        std::getline(in, s);
        ASSERT_EQ(s, correctFile[i]);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
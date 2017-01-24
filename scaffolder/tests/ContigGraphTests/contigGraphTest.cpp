#include "ContigGraph/ContigGraph.h"
#include "gtest/gtest.h"

class ContigGraphTest : public  ::testing::Test {
};

TEST_F(ContigGraphTest, addNewLib) {
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
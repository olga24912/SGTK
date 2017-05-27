#ifndef SCAFFOLDER_BLOCKSPLITDOTWRITER_H
#define SCAFFOLDER_BLOCKSPLITDOTWRITER_H

#include "DotWriter.h"

class BlockSplitDotWriter : public DotWriter {
public:
    BlockSplitDotWriter(Filter *filter, FileValidator *validator, int maxVert, int maxEdge);
protected:
    const int minBlockDist = 100;

    struct vertBlock {
        int vertId;
        int coordB;
        int coordE;
        int hasOutsideEdge = 0;
        vertBlock(){}
        vertBlock(int vertId, int coordB, int coordE) : vertId(vertId), coordB(coordB), coordE(coordE) {}

        bool operator < (vertBlock b2) {
            return coordB > b2.coordB;
        }
    };

    void writeOneVertexSet(std::vector<int> vert, std::string fileName) override;

    std::vector<vertBlock> splitOnBlocks(int i);

    void findOutsideEdge(std::vector<std::vector<vertBlock>>& blocks);

    bool isOutsideEdge(int e, const std::vector<std::vector<vertBlock>> &vector);

    void writeOneVertBlock(std::vector<vertBlock> &bl, std::ofstream &out);

    void writeBlockEdges(const std::vector<vertBlock> &bl, std::ofstream &out);

    void writeEdges(const std::vector<std::vector<vertBlock>> &bl, std::ofstream &out);

    int findBlockId(int v, int e, const std::vector<std::vector<vertBlock>> &bl);
};


#endif //SCAFFOLDER_BLOCKSPLITDOTWRITER_H

#ifndef SCAFFOLDER_SCAFFOLDS_H
#define SCAFFOLDER_SCAFFOLDS_H

#include <string>
#include <vector>

class Scaffolds {
private:
    class Node {
    public:
        Node* next;
        Node* priv;
        int id;
        Node* rev;
    };

    const int GAP_SIZE = 100;

    std::vector<std::string> contigs;
    std::vector<std::string> contigsName;
    std::vector<Node> contigsNode;
    std::vector<Node*> scaffolds;

    void saveContigs(std::string contigFile);
    std::string createRevCompl(std::string s);

public:
    Scaffolds(std::string contigFile);
    void print(std::string out);
    void addConnection(int id1, int id2);
    void brokeConnection(int id1);
    void brokeConnectionTo(int id2);

    int lineId(int id);
};


#endif //SCAFFOLDER_SCAFFOLDS_H

#include "RNAPairReadGraphBuilder.h"

std::string RNAPairReadGraphBuilder::getLibColor() {
    int color[3] = {rand()%100, rand()%150, 255};
    return colorToString(color);
}

//
// Created by olga on 08.10.16.
//

#include "RNAPairReadGraphBuilder.h"

string RNAPairReadGraphBuilder::getLibColor() {
    int color[3] = {rand()%100, rand()%150, 255};
    return GraphUtils::colorToString(color);
}

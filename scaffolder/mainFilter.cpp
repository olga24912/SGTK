#include <iostream>
#include "Builder/GraphControl.h"

using namespace std;

#include "Filter/InteractiveFilter.h"

int main(int argc, char **argv) {
    InteractiveFilter filter;
    filter.main();
    return 0;
}
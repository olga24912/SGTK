include_directories(.)

set(SOURCE_FILES ExonInfo.cpp ExonInfo.h ExonsInfo.cpp ExonsInfo.h)

add_library(findExon_lib ${SOURCE_FILES})
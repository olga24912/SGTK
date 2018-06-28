find_package (ZLIB)
find_package (BZip2)

include_directories(.)

set(SOURCE_FILES ReadsSplitter.cpp ReadsSplitter.h
        ReadsSplitter50.cpp ReadsSplitter50.h
        SplitterByUnmappedEnd.cpp SplitterByUnmappedEnd.h
        Utils/SeqanUtils.cpp Utils/SeqanUtils.h)

add_library(read_splitter_lib ${SOURCE_FILES})
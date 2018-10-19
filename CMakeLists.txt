cmake_minimum_required(VERSION 3.5)
project(SGTK)

# Define option for static / dynamic build.
option(SGTK_STATIC_BUILD "Link SGTK statically" OFF)
if (SGTK_STATIC_BUILD)
    # it'll make cmake to find libraries archives, not dynamic link
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    set(LINK_SEARCH_START_STATIC TRUE)
    set(LINK_SEARCH_END_STATIC TRUE)
    # This is dirty hack to get rid of -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS "-Wl,-Bstatic")
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS "-Wl,-Bstatic")

    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
        add_definitions(-static)
    endif()

    set(Boost_USE_STATIC_LIBS        ON)
    set(Boost_USE_STATIC_RUNTIME     ON)
endif()


find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -std=c++14 -O3 -DSEQAN_HAS_ZLIB=1")

include_directories(src)

add_subdirectory(src)
#add_subdirectory(tests)

# Main pipeline script
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/rna_scaffolder.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/sgtk.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/ref.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/contigs.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/read_1.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/read_2.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/test_dataset/scaf.info"
        DESTINATION share/test_dataset
        COMPONENT runtime)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape-dagre.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/dagre.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape-qtip.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js-navigator/cytoscape-navigator.js"
        DESTINATION bin/scripts/external/cytoscape.js-navigator/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js-navigator/bower.json"
        DESTINATION bin/scripts/external/cytoscape.js-navigator/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js-navigator/gulpfile.js"
        DESTINATION bin/scripts/external/cytoscape.js-navigator/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js-navigator/package.json"
        DESTINATION bin/scripts/external/cytoscape.js-navigator/
        COMPONENT runtime)


install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/scaffoldgraph.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/showListScript.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/findAreaScript.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/script.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/search.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleAlongChromosomes.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleScaffolds.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleDiffInLibs.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/defaultWeight.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/zoomChoose.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleAmbiguous.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/freeLayout.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/chromFindVertexRank.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/Icon.png"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/light.css"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/mainPage.html"
        DESTINATION bin/
        COMPONENT runtime)

cmake_minimum_required(VERSION 3.3)
project(scaffolder)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -std=c++14 -O3 -DSEQAN_HAS_ZLIB=1")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

# Main pipeline script
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/rna_scaffolder.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/visualize_scaffold_graph.py"
        DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/ref.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/contigs.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/read_1.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/read_2.fasta"
        DESTINATION share/test_dataset
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/resources/scaf.info"
        DESTINATION share/test_dataset
        COMPONENT runtime)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape-dagre.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/vivagraph.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/dagre.js"
        DESTINATION bin/scripts/external/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/external/cytoscape-qtip.js"
        DESTINATION bin/scripts/external/
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
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleAlongChromosomes.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/habdleScaffolds.js"
        DESTINATION bin/scripts/
        COMPONENT runtime)
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/UI/scripts/handleDiffInLibs.js"
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

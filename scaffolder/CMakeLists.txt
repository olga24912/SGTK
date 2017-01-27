cmake_minimum_required(VERSION 3.3)
project(scaffolder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(build mainBuildGraph.cpp tests/BuilderTests/builderTest.cpp)
add_executable(filter mainFilter.cpp)

target_link_libraries(build scaffolder_lib)
target_link_libraries(filter scaffolder_lib)
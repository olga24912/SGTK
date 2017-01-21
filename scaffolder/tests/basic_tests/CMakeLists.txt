include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runBasicScaffolderTests
        basic_check.cpp)

target_link_libraries(runBasicScaffolderTests gtest gtest_main)
target_link_libraries(runBasicScaffolderTests scaffolder_lib)
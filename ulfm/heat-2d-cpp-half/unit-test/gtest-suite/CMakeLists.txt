cmake_minimum_required(VERSION 3.9.1)
 
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(test-suite test-suite.cpp)

target_link_libraries(test-suite ${GTEST_LIBRARIES} pthread)

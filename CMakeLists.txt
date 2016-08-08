cmake_minimum_required(VERSION 3.5)
project(Elasticity)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#Find and setup deal ii
find_package (deal.II REQUIRED
              HINTS $ENV{DEAL_II_DIR})
deal_ii_initialize_cached_variables()

file (GLOB_RECURSE LIBRARY_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/libelast/*.cpp")
add_library (eLast SHARED ${LIBRARY_SOURCES})
deal_ii_setup_target (eLast)

file (GLOB_RECURSE LAUNCHER_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/*.cpp")
add_library (launcher ${LIBRARY_SOURCES})
deal_ii_setup_target (launcher)

set(SOURCE_FILES main.cpp)
add_executable(Elasticity ${SOURCE_FILES})
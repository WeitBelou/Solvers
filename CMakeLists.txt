cmake_minimum_required(VERSION 3.5)
project(Solvers)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set (CMAKE_INCLUDE_CURRENT_DIR ON)

#Boost
find_package(Boost COMPONENTS
             program_options
             filesystem
             REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

#Find and setup deal ii
find_package (deal.II REQUIRED
              HINTS $ENV{DEAL_II_DIR})
deal_ii_initialize_cached_variables()

file (GLOB_RECURSE LIBRARY_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/libsolvers/*.cpp")
add_library (solvers SHARED ${LIBRARY_SOURCES})
deal_ii_setup_target (solvers)

file (GLOB_RECURSE LAUNCHER_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/*.cpp")
add_library (launcher ${LAUNCHER_SOURCES})
target_link_libraries (launcher ${Boost_FILESYSTEM_LIBRARIES})
deal_ii_setup_target (launcher)

set(SOURCE_FILES main.cpp)
add_executable(Solvers ${SOURCE_FILES})
target_link_libraries (Solvers ${Boost_LIBRARIES})
target_link_libraries(Solvers solvers launcher)

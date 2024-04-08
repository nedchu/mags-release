cmake_minimum_required(VERSION 2.8)

project(core-sum)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -O3 -fopenmp")
find_package(OpenMP REQUIRED)

include_directories(./src)
aux_source_directory(./src INCLUDES)
add_library(INCLUDES_LIB ${INCLUDES})

add_executable(mags ./run/run_mags.cpp)
target_link_libraries(mags INCLUDES_LIB)
add_executable(mags_dm ./run/run_mags_dm.cpp)
target_link_libraries(mags_dm INCLUDES_LIB)
add_executable(pmags_dm ./run/run_para_mags_dm.cpp)
target_link_libraries(pmags_dm INCLUDES_LIB)
add_executable(pmags ./run/run_para_mags.cpp)
target_link_libraries(pmags INCLUDES_LIB)

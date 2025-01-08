cmake_minimum_required(VERSION 3.16)

project(sorting_exercises)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(include)
add_executable(Sort source/sorting_timing.cpp)
#target_compile_options(Sort PRIVATE -O1)
cmake_minimum_required(VERSION 3.16)

project(matrix_exercises)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(include)
add_executable(MatrixOps source/matrix_operations.cpp)
# target_compile_options(MatrixOps PRIVATE -O3)
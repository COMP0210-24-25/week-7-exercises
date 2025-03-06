
add_executable(Sum Summation.cpp)
target_include_directories(Sum PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_options(Sum PRIVATE -mavx)  # remove -mavx if your cpu doesn't support it
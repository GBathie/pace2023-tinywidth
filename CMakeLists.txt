cmake_minimum_required(VERSION 3.15)

project(main)

file(GLOB SRC src/*.cpp)
add_executable(main ${SRC})
add_compile_definitions(OPTIL)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -mavx2 -std=c++2a")
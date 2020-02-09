cmake_minimum_required(VERSION 2.8.12)
project(Homework5 C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fopenmp -lm -std=c11")

set(CMAKE_C_COMPILER gcc)

#link_directories( /usr/local/opt/libomp/lib )
#include_directories( /usr/local/opt/libomp/include/ )

add_executable(Homework5-parallel parallel.c)
#target_link_libraries(Homework5-parallel omp)

add_executable(Homework5-serial main.c)

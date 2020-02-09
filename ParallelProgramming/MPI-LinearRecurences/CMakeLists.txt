cmake_minimum_required(VERSION 2.8.12)
project(Homework3 C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -lm")

set(CMAKE_C_COMPILER mpicc)

add_executable(Homework3-CR main.c)
add_executable(Homework3-RD rd.c)


cmake_minimum_required(VERSION 2.8.12)
project(Homework2 C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -lm -D FOX ")

set(CMAKE_C_COMPILER mpicc)

add_executable(Homework2 main.c)


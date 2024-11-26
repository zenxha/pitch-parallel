cmake_minimum_required(VERSION 3.10)
project(pitch-parallel)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(pitch-parallel src/main.cpp)
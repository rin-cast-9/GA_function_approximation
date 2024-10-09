cmake_minimum_required(VERSION 3.16)

# project name and version
project(genetic_algorithm VERSION 1.0 LANGUAGES CXX)

# c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# source files
add_executable(main ./src/main.cpp)
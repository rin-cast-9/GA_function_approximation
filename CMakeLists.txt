cmake_minimum_required(VERSION 3.16)

# project name and version
project(genetic_algorithm VERSION 1.0 LANGUAGES CXX)

# c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_definitions(-DPROJECT_ROOT_DIR="${CMAKE_SOURCE_DIR}")

# boost installation path
set(BOOST_ROOT "/opt/homebrew/Cellar/boost/1.86.0")
set(Boost_NO_SYSTEM_PATHS TRUE)

# find boost libraries
find_package(Boost 1.86 REQUIRED COMPONENTS filesystem)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(main ./src/main.cpp)
    target_link_libraries(main ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found!")
endif()
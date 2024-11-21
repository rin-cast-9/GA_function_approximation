cmake_minimum_required(VERSION 3.16)

# project name and version
project(genetic_algorithm VERSION 1.0 LANGUAGES CXX)

# c++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_definitions(-DPROJECT_ROOT_DIR="${CMAKE_SOURCE_DIR}")
add_definitions(-DPROJECT_BINARY_DIR="${CMAKE_BINARY_DIR}")

# boost installation path
set(BOOST_ROOT "/opt/homebrew/Cellar/boost/1.86.0")
set(Boost_NO_SYSTEM_PATHS TRUE)

# find boost libraries
find_package(Boost 1.86 REQUIRED COMPONENTS filesystem)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    file (GLOB SOURCES "./src/*.cpp")
    file (GLOB HEADERS "./src/*.hpp")

    # define the custom command for compiling the metal file
    add_custom_command(
        OUTPUT ${CMAKE_BINARY_DIR}/default.metallib
        COMMAND xcrun -sdk macosx metal -c ${CMAKE_SOURCE_DIR}/src/crossover.metal -o ${CMAKE_BINARY_DIR}/crossover.air
        COMMAND xcrun -sdk macosx metal -c ${CMAKE_SOURCE_DIR}/src/fitness.metal -o ${CMAKE_BINARY_DIR}/fitness.air
        COMMAND xcrun -sdk macosx metal -c ${CMAKE_SOURCE_DIR}/src/mutation.metal -o ${CMAKE_BINARY_DIR}/mutation.air
        COMMAND xcrun -sdk macosx metallib ${CMAKE_BINARY_DIR}/crossover.air ${CMAKE_BINARY_DIR}/fitness.air ${CMAKE_BINARY_DIR}/mutation.air -o ${CMAKE_BINARY_DIR}/default.metallib
        DEPENDS ${CMAKE_SOURCE_DIR}/src/crossover.metal ${CMAKE_SOURCE_DIR}/src/fitness.metal
        COMMENT "Compiling metal shader crossover.metal and fitness.metal into default.metallib"
        VERBATIM
    )

    add_custom_target(
        metal_lib ALL
        DEPENDS ${CMAKE_BINARY_DIR}/default.metallib
    )

    add_executable(main ${SOURCES})
    add_dependencies(main metal_lib)
    target_include_directories(main PRIVATE /usr/local/include/metal-cpp)
    target_link_libraries(main ${Boost_LIBRARIES})
    target_link_libraries(main "-framework Metal" "-framework Foundation" "-framework QuartzCore")
else()
    message(FATAL_ERROR "Boost not found!")
endif()

cmake_minimum_required(VERSION 3.5)

project(integration_algorithms)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(!COMPILER_SUPPORTS_CXX17)
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} does not support C++17. Please use a different C++ compiler.")
    return()
endif()

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++2a -mavx2 -mfma -O3 -Wno-sign-conversion")
SET(GCC_COVERAGE_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

# Create a sources variable with a link to all cpp files to compile
set(integration_algorithms_SOURCES
    integration_algorithms.cpp
)

add_executable(integration_algorithms ${integration_algorithms_SOURCES})

target_include_directories(integration_algorithms
    PRIVATE 
)

cmake_minimum_required(VERSION 3.5)

project(test_integration)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(!COMPILER_SUPPORTS_CXX17)
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} does not support C++17. Please use a different C++ compiler.")
    return()
endif()

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++17 -mavx2 -mfma -O3 -Wno-sign-conversion")
SET(GCC_COVERAGE_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

# Create a sources variable with a link to all cpp files to compile
set(test_integration_SOURCES
    test_integration.cpp
    libs/numerical_algorithms.cpp
)

add_executable(test_integration ${test_integration_SOURCES})

target_include_directories(test_integration
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

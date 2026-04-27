cmake_minimum_required(VERSION 3.16)
project(linal_lab1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(linal_core
    matrix.cpp
    gauss_default.cpp
    gauss_advanced.cpp
    LU_decomposition.cpp
)

target_include_directories(linal_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (MSVC)
    target_compile_options(linal_core PRIVATE /W4)
else()
    target_compile_options(linal_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(linal_lab1
    tests.cpp
)

target_link_libraries(linal_lab1 PRIVATE linal_core)

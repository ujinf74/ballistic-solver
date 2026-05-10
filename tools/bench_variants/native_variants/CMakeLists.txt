cmake_minimum_required(VERSION 3.20)

project(bench_variants_core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(bench_variants_core MODULE
    bench_variants_module.cpp
)

target_include_directories(bench_variants_core PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/../../../include
)

set_target_properties(bench_variants_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

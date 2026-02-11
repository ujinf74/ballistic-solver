cmake_minimum_required(VERSION 3.20)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

set(BS_VERSION "${SKBUILD_PROJECT_VERSION}")

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../.. ${CMAKE_BINARY_DIR}/bs_root)

pybind11_add_module(_core MODULE
  pybind_module.cpp
)

target_link_libraries(_core PRIVATE ballistic_solver)

target_include_directories(_core PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/../../src
)

set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ballistic_solver"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ballistic_solver"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ballistic_solver"
)

install(TARGETS _core
  LIBRARY DESTINATION ballistic_solver
  RUNTIME DESTINATION ballistic_solver
  ARCHIVE DESTINATION ballistic_solver
)

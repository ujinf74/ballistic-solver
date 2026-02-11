cmake_minimum_required(VERSION 3.20)
if(DEFINED SKBUILD_PROJECT_VERSION AND DEFINED SKBUILD_PROJECT_NAME)
  project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)
else()
  project(ballistic_solver_py LANGUAGES C CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE
    pybind_module.cpp
)

if(DEFINED SKBUILD_PROJECT_VERSION)
  target_compile_definitions(_core PRIVATE
    BALLISTIC_SOLVER_VERSION_STRING="${SKBUILD_PROJECT_VERSION}"
  )
endif()

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

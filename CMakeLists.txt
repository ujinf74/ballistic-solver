cmake_minimum_required(VERSION 3.20)

if(DEFINED BS_VERSION)
  set(BS_VERSION_DETECTED "${BS_VERSION}")
else()
  bs_detect_version_from_git(BS_VERSION_DETECTED)
endif()

project(ballistic_solver VERSION "${BS_VERSION_DETECTED}" LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(ballistic_solver SHARED
    src/ballistic_solver_c_api.cpp
)

target_include_directories(ballistic_solver
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(ballistic_solver PRIVATE
    BALLISTIC_SOLVER_EXPORTS
    BALLISTIC_SOLVER_ABI_VERSION=2
    BALLISTIC_SOLVER_VERSION_STRING=\"${PROJECT_VERSION}\"
)

if (MSVC)
    target_compile_options(ballistic_solver PRIVATE /W4 /utf-8)
else()
    target_compile_options(ballistic_solver PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(ballistic_solver_example
    examples/c/basic.c
)
target_link_libraries(ballistic_solver_example PRIVATE ballistic_solver)

include(CTest)
if (BUILD_TESTING)
    add_executable(smoke_test
        tests/smoke_test.c
    )
    target_link_libraries(smoke_test PRIVATE ballistic_solver)

    add_test(
        NAME smoke_test
        COMMAND $<TARGET_FILE:smoke_test>
    )
endif()

include(GNUInstallDirs)

install(TARGETS ballistic_solver
    EXPORT ballistic_solverTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Minimal find_package(CONFIG) support
set(_pkgdir ${CMAKE_INSTALL_LIBDIR}/cmake/ballistic_solver)

install(EXPORT ballistic_solverTargets
    FILE ballistic_solverTargets.cmake
    NAMESPACE ballistic_solver::
    DESTINATION ${_pkgdir}
)

include(CMakePackageConfigHelpers)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/ballistic_solverConfig.cmake"
"include(\"\${CMAKE_CURRENT_LIST_DIR}/ballistic_solverTargets.cmake\")\n"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ballistic_solverConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ballistic_solverConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ballistic_solverConfigVersion.cmake"
    DESTINATION ${_pkgdir}
)

set_target_properties(ballistic_solver PROPERTIES OUTPUT_NAME "ballistic_solver")

if(DEFINED SKBUILD_PLATLIB_DIR)
  set(BS_NATIVE_DIR "${SKBUILD_PLATLIB_DIR}/ballistic_solver/_native")

  install(TARGETS ballistic_solver
    RUNTIME DESTINATION "${BS_NATIVE_DIR}"  # Windows .dll
    LIBRARY DESTINATION "${BS_NATIVE_DIR}"  # Linux .so / macOS .dylib
    ARCHIVE DESTINATION "${BS_NATIVE_DIR}"  # .lib
  )
endif()


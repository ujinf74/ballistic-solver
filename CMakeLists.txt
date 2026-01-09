cmake_minimum_required(VERSION 3.16)

project(ballistic_solver
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(ballistic INTERFACE)
target_include_directories(ballistic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (MSVC)
    target_compile_options(ballistic INTERFACE /W4 /utf-8)
else()
    target_compile_options(ballistic INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_executable(ballistic_example examples/basic.cpp)
target_link_libraries(ballistic_example PRIVATE ballistic)

add_library(ballistic_c SHARED
    src/ballistic_c_api.cpp
)

target_link_libraries(ballistic_c
    PRIVATE
        ballistic
)

if (MSVC)
    target_compile_options(ballistic_c PRIVATE /W4 /utf-8)
else()
    target_compile_options(ballistic_c PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if (BUILD_TESTING)
    add_executable(ballistic_smoke_test tests/smoke_test.cpp)
    target_link_libraries(ballistic_smoke_test PRIVATE ballistic)
    add_test(
        NAME ballistic_smoke
        COMMAND $<TARGET_FILE:ballistic_smoke_test>
    )
endif()
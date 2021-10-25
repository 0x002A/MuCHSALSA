cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(VERSION_MAJOR 0 CACHE STRING "Project major version number.")
set(VERSION_MINOR 0 CACHE STRING "Project minor version number.")
set(VERSION_PATCH 1 CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

set(CMAKE_C_COMPILER "clang-12")
set(CMAKE_CXX_COMPILER "clang++-12")

set(TARGET "muchsalsa")
set(LIBNAME "ms")

project(${TARGET} CXX)
include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)

add_library(project_options INTERFACE)

option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
if (ENABLE_BUILD_WITH_TIME_TRACE)
    target_compile_options(project_options INTERFACE -ftime-trace)
endif ()

add_library(project_warnings INTERFACE)

include(cmake/Cache.cmake)

include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

include(cmake/StaticAnalyzers.cmake)

option(ENABLE_TESTING "Enable Test Builds" ON)

option(ENABLE_PCH "Enable Precompiled Headers" ON)
if (ENABLE_PCH)
    target_precompile_headers(
            project_options
            INTERFACE
            <cstddef>
            <map>
            <string>
            <tuple>
            <unordered_map>
            <utility>
            <vector>)
endif ()

if (ENABLE_TESTING)
    enable_testing()
endif ()

add_subdirectory("lib${LIBNAME}")

file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")
add_executable(${TARGET} ${SOURCES})

option(ENABLE_TRACK_MEMORY_USAGE "Enable tracking of memory allocations" OFF)
if (ENABLE_TRACK_MEMORY_USAGE)
    target_compile_definitions(${TARGET} PRIVATE TRACK_MEMORY_USAGE)
endif ()

option(ENABLE_BUILD_FAT_EXE "Enable linking dependencies statically" OFF)
if (ENABLE_BUILD_FAT_EXE)
    target_link_options(${TARGET} PRIVATE -static)
endif ()

target_include_directories(${TARGET} PRIVATE "include")
target_link_libraries(
        ${TARGET}
        PRIVATE ${LIBNAME}
        PRIVATE c++
        PRIVATE c++abi
        PRIVATE project_options
        PRIVATE project_warnings)

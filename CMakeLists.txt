cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_C_COMPILER "clang-10")
set(CMAKE_CXX_COMPILER "clang++-10")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

set(TARGET "lazybastard")
set(LIBNAME "lb")

project(${TARGET} CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=*,-fuchsia-default-arguments -extra-arg=-std=c++17)

enable_testing()

add_subdirectory("lib${LIBNAME}")

file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" "src/*.h")

add_executable(${TARGET} ${SOURCES})
target_include_directories(${TARGET} PUBLIC "include")
target_link_libraries(${TARGET} ${LIBNAME})

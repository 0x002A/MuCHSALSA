cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)

configure_file(${doxyfile_in} ${doxyfile} @ONLY)

set(TARGET "lb")

file(GLOB_RECURSE SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/lib${TARGET}" "src/*.cpp" "src/*.h")

add_library(${TARGET}  STATIC  ${SOURCES})
target_include_directories(${TARGET} PUBLIC "${CMAKE_SOURCE_DIR}/include/${TARGET}")

add_custom_target(doc
    COMMAND doxygen ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

project("${TARGET}_test" CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE TEST_SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/lib${TARGET}" "tests/*.cpp" "tests/*.h")

add_executable("${TARGET}_test" ${TEST_SOURCES})
add_test(NAME "${TARGET}_test" COMMAND "${TARGET}_test")

target_include_directories("${TARGET}_test" PUBLIC "${CMAKE_SOURCE_DIR}/include/${TARGET}")
target_link_libraries("${TARGET}_test" ${TARGET} gtest)

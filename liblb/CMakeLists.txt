cmake_minimum_required(VERSION 3.17.0)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include(FetchContent)

FetchContent_Declare(
  gsl
  GIT_REPOSITORY https://github.com/microsoft/GSL.git
  GIT_TAG        v3.0.1
)

FetchContent_GetProperties(gsl)
if(NOT gsl_POPULATED)
  FetchContent_Populate(gsl)
endif()

include(../cmake/Doxygen.cmake)

set(TARGET "lb")

file(GLOB_RECURSE SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/lib${TARGET}" "src/*.cpp")

add_library(${TARGET} STATIC ${SOURCES})

target_include_directories(${TARGET} PUBLIC "${CMAKE_SOURCE_DIR}/include/${TARGET}")
target_include_directories(${TARGET} PUBLIC "${gsl_SOURCE_DIR}/include")

target_link_libraries(${TARGET}
        PUBLIC Threads::Threads
        PRIVATE project_options
        project_warnings)

set(TEST_SRC_DIR "${CMAKE_SOURCE_DIR}/lib${TARGET}")
include(../cmake/Testing.cmake)
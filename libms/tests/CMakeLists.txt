FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

FetchContent_GetProperties(googletest)
if (NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

set_target_properties(
        "gtest"
        PROPERTIES
        CXX_CLANG_TIDY ""
)
set_target_properties(
        "gtest_main"
        PROPERTIES
        CXX_CLANG_TIDY ""
)
set_target_properties(
        "gmock"
        PROPERTIES
        CXX_CLANG_TIDY ""
)
set_target_properties(
        "gmock_main"
        PROPERTIES
        CXX_CLANG_TIDY ""
)

file(GLOB_RECURSE TEST_SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/lib${TARGET}/tests" "*.cpp")

add_executable("${TARGET}_test" ${TEST_SOURCES})
set_target_properties(
        "${TARGET}_test"
        PROPERTIES
        CXX_CLANG_TIDY ""
)

target_link_libraries("${TARGET}_test" PRIVATE project_options gtest ${TARGET})

add_test(NAME "${TARGET}_test" COMMAND "${TARGET}_test" "${CMAKE_SOURCE_DIR}/test_data")

option(ENABLE_TESTING "Enable unit-tests" ON)

if (ENABLE_TESTING)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.10.0
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

    project("${TARGET}_test" CXX)

    file(GLOB_RECURSE TEST_SOURCES RELATIVE ${TEST_SRC_DIR} "tests/*.cpp")

    add_executable("${TARGET}_test" ${TEST_SOURCES})
    set_target_properties(
            "${TARGET}_test"
            PROPERTIES
            CXX_CLANG_TIDY ""
    )

    target_link_libraries("${TARGET}_test" ${TARGET} gtest)

    enable_testing()
    add_test(NAME "${TARGET}_test" COMMAND "${TARGET}_test")
endif ()
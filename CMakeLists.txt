cmake_minimum_required(VERSION 3.31)
project(STABAAL LANGUAGES CXX)

if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    include(CheckIPOSupported)
    check_ipo_supported(LANGUAGES CXX)
endif()

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(stabaal STATIC ${SOURCES})
target_compile_features(stabaal PUBLIC cxx_std_20)

target_include_directories(stabaal PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

option(BUILD_BENCHMARKS "Build the tableau benchmark" OFF)
if(BUILD_BENCHMARKS)
    add_executable(stabaal_bench benchmarks/bench.cpp)
    target_link_libraries(stabaal_bench PRIVATE stabaal)
endif()

include(CTest)
if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.1
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(stabaal_tests tests/test_tableau.cpp tests/test_dynamic_bitset.cpp)
    target_link_libraries(stabaal_tests PRIVATE stabaal Catch2::Catch2WithMain)

    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)
    catch_discover_tests(stabaal_tests)
endif()

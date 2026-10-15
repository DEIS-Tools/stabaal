cmake_minimum_required(VERSION 3.16)
project(STABAAL LANGUAGES CXX)

find_package(Boost REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(stabaal_lib
    src/tableau.cpp
)
target_include_directories(stabaal_lib PUBLIC include)

add_executable(stabaal main.cpp)
target_link_libraries(stabaal PRIVATE stabaal_lib Boost::headers)

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(Catch2)

enable_testing()
add_executable(stabaal_tests tests/test_tableau.cpp)
target_link_libraries(stabaal_tests PRIVATE stabaal_lib Boost::headers Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(stabaal_tests)

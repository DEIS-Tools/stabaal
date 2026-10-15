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

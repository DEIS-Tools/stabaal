cmake_minimum_required(VERSION 3.28)
project(STABAAL LANGUAGES CXX)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

include(FetchContent)

if(CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    include(CheckIPOSupported)
    check_ipo_supported(LANGUAGES CXX)
endif()

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX "/src/bindings/")
add_library(stabaal STATIC ${SOURCES})
target_compile_features(stabaal PUBLIC cxx_std_20)

target_include_directories(stabaal PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

option(BUILD_BENCHMARKS "Build the benchmark" OFF)
if(BUILD_BENCHMARKS)
    add_executable(stabaal_bench benchmarks/bench.cpp)
    target_link_libraries(stabaal_bench PRIVATE stabaal)
endif()

option(BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
if(BUILD_PYTHON_BINDINGS)
    set_target_properties(stabaal PROPERTIES POSITION_INDEPENDENT_CODE ON)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    find_package(pybind11 3.1 CONFIG REQUIRED)

    pybind11_add_module(stabaal_python src/bindings/python/bindings.cpp)
    target_link_libraries(stabaal_python PRIVATE stabaal)
    set_target_properties(stabaal_python PROPERTIES OUTPUT_NAME stabaal)
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/stubs/stabaal.pyi"
        COMMAND "${Python_EXECUTABLE}" -m pybind11_stubgen stabaal
                -o "${CMAKE_CURRENT_BINARY_DIR}/stubs" --exit-code
        WORKING_DIRECTORY "$<TARGET_FILE_DIR:stabaal_python>"
        DEPENDS stabaal_python
        VERBATIM
    )
    add_custom_target(stabaal_stubs ALL
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/stubs/stabaal.pyi")
    install(TARGETS stabaal_python LIBRARY DESTINATION . RUNTIME DESTINATION .)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/stubs/stabaal.pyi" DESTINATION .)
endif()

include(CTest)
if(BUILD_TESTING)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.1
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(stabaal_tests tests/test_inverse_tableau.cpp tests/test_pauli_string.cpp)
    target_link_libraries(stabaal_tests PRIVATE stabaal Catch2::Catch2WithMain)

    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)
    catch_discover_tests(stabaal_tests)
endif()

foreach(target IN ITEMS stabaal stabaal_bench stabaal_python stabaal_tests)
    if(TARGET ${target})
        target_compile_options(${target} PRIVATE
            "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra>")
    endif()
endforeach()

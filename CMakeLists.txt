cmake_minimum_required(VERSION 3.14)
project(RSA)

set(CMAKE_CXX_STANDARD 17)
# Compiler Flags
set(GCC_WARNINGS "-Wall -Wextra -Wpedantic -Wunused -Wformat=2 -Wswitch-enum -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wcast-align -Wcast-qual -Wconversion -Wsign-conversion -Wzero-as-null-pointer-constant -Wlogical-op -Wmissing-declarations -Wctor-dtor-privacy -Wnoexcept -Wnon-virtual-dtor -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wfloat-equal -Wredundant-decls -Wshadow -Wwrite-strings -Wno-long-long")
set(CLANG_WARNINGS "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-documentation-unknown-command")


set(CLANG_DEBUG "-O0 -g")
set(CLANG_RELEASE "-O3")

set(GCC_DEBUG "-O0 -g")
set(GCC_RELEASE "-O3")


# Set up the compiler flags used
if(CMAKE_COMPILER_IS_GNUCXX)
    # gcc
    set(CMAKE_CXX_FLAGS "${GCC_WARNINGS}")

    set(CMAKE_CXX_FLAGS_DEBUG "${GCC_DEBUG}")
    set(CMAKE_CXX_FLAGS_RELEASE "${GCC_RELEASE}")
else()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        # clang
        set(CMAKE_CXX_FLAGS "${CLANG_WARNINGS}")

        set(CMAKE_CXX_FLAGS_DEBUG "${CLANG_DEBUG}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CLANG_RELEASE}")
    else()
        # Generic flags for undetected compiler
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    endif()
endif()

set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(GMP REQUIRED)

#if (!GMP_FOUND)
#    message(FATAL_ERROR "GMP not found.")
#endif()

#include_directories(${GMP_INCLUDE_DIR})

# This is for Mac Homebrew version of GMP. Likely needs to change for other setups
#link_directories("/usr/local/opt/gmp/lib/")

#message(FATAL_ERROR "Found at ${GMP_LIBRARY}")
#message(FATAL_ERROR "Found at ${GMPXX_LIBRARY}")

add_executable(RSA RSA.cpp)
target_link_libraries(RSA GMP::GMPXX GMP::GMP)


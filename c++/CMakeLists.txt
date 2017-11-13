cmake_minimum_required(VERSION 2.8.12)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/msvc-specific.cmake)

SET(PROJECT_NAME BoundednessMPCSolver)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "https://gite.lirmm.fr/adrien_escande_at_gmail_com/BoundednessMPCSolver")
SET(PROJECT_DEBUG_POSTFIX "_d")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

SETUP_PROJECT()

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-conversion -std=c++0x -pedantic")
endif()

option(PYTHON_BINDINGS "Generate Python bindings" ON)
option(STATISTICS "Use statistics" ON)

#########################
# External dependencies
SEARCH_FOR_EIGEN()

SET(BOOST_COMPONENTS python system)
IF(PYTHON_BINDINGS)
    set(BOOST_COMPONENTS python unit_test_framework timer system)
ELSE()
    set(BOOST_COMPONENTS unit_test_framework timer system)
ENDIF()

SEARCH_FOR_BOOST()

#########################

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
IF(MSVC)
  SET(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
ENDIF(MSVC)

add_subdirectory(src)
add_subdirectory(tests)

# Add dependency towards the library in the pkg-config file.
PKG_CONFIG_APPEND_LIBS(BoundednessMPCSolver)

IF(PYTHON_BINDINGS)
    # Add dependency towards the python bindings in the pkg-config file.
    add_subdirectory(bindings)
    PKG_CONFIG_APPEND_LIBS(PyBoundednessMPCSolver)
ENDIF()

SETUP_PROJECT_FINALIZE()

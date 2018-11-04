# Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
#
# This file is part of CPS.
#
# CPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CPS.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 2.8.12)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/msvc-specific.cmake)

SET(PROJECT_NAME CaptureProblemSolver)
SET(PROJECT_DESCRIPTION "Solver for capture problems when walking with variable height")
SET(PROJECT_URL "https://github.com/jrl-umi3218/CaptureProblemSolver")
SET(PROJECT_DEBUG_POSTFIX "_d")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

SETUP_PROJECT()

if(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-conversion -std=c++0x -pedantic")
endif()

option(PYTHON_BINDINGS "Generate Python bindings" ON)
option(STATISTICS "Use statistics" OFF)

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

if(STATISTICS)
  add_definitions(-DUSE_STATS)
endif()

add_subdirectory(src)
add_subdirectory(tests)

# Add dependency towards the library in the pkg-config file.
PKG_CONFIG_APPEND_LIBS(CaptureProblemSolver)

IF(PYTHON_BINDINGS)
    # Add dependency towards the python bindings in the pkg-config file.
    add_subdirectory(bindings)
ENDIF()

SETUP_PROJECT_FINALIZE()

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

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${Boost_INCLUDE_DIRS})
if(WIN32)
  #This is one of the way to avoid link errors related to static variables in program_options
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif(WIN32)

add_definitions(-DTESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

macro(addUnitTest name)
  add_executable(${name} ${name}.cpp ${ARGN})
  target_link_libraries(${name} ${Boost_LIBRARIES} CaptureProblemSolver)
  add_test(${name}Unit ${name})
  # Adding a project configuration file (for MSVC only)
  GENERATE_MSVC_DOT_USER_FILE(${name})
endmacro(addUnitTest)

addUnitTest(MatrixComputationsTest)
addUnitTest(LeastSquareTest)
addUnitTest(LinearConstraintsTest)
addUnitTest(ProblemTest)
addUnitTest(ObjectiveTest)
addUnitTest(SQPTest SQPTestCommon.h SQPTestCommon.cpp)

set(TESTSOURCES
  main.cpp
  SQPTestCommon.cpp
  timings.cpp
)

set(TESTHEADERS
  SQPTestCommon.h
  timings.h
)

add_executable(main ${TESTSOURCES} ${TESTHEADERS})
target_link_libraries(main PUBLIC CaptureProblemSolver)
GENERATE_MSVC_DOT_USER_FILE(main)


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

INCLUDE(../cmake/python.cmake)
FINDPYTHON()
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

set(MODULE_NAME "CaptureProblemSolver")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}/__init__.py COPYONLY)

set(PYCPS_SOURCES PyCaptureProblemSolver.cpp)

set(CPS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(PYCPS_HEADERS converters.h ${CPS_INCLUDE_DIR}/cps/cps_api.h)

add_library(PyCaptureProblemSolver SHARED ${PYCPS_SOURCES} ${PYCPS_HEADERS})
target_link_libraries(PyCaptureProblemSolver CaptureProblemSolver boost_numpy)
TARGET_LINK_BOOST_PYTHON(PyCaptureProblemSolver)
set_target_properties(PyCaptureProblemSolver PROPERTIES PREFIX "")
install(TARGETS PyCaptureProblemSolver DESTINATION "${PYTHON_SITELIB}/${MODULE_NAME}")

PYTHON_INSTALL_BUILD(${MODULE_NAME} __init__.py "${PYTHON_SITELIB}")

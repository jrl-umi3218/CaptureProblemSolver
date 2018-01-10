INCLUDE(../cmake/python.cmake)
FINDPYTHON()
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

set(MODULE_NAME "BoundednessMPCSolver")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}/__init__.py COPYONLY)

set(PYBMS_SOURCES PyBoundednessMPCSolver.cpp)

set(BMS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(PYBMS_HEADERS converters.h ${BMS_INCLUDE_DIR}/bms/bms_api.h)

add_library(PyBoundednessMPCSolver SHARED ${PYBMS_SOURCES} ${PYBMS_HEADERS})

target_link_libraries(PyBoundednessMPCSolver BoundednessMPCSolver boost_numpy)
TARGET_LINK_BOOST_PYTHON(PyBoundednessMPCSolver)

set_target_properties(PyBoundednessMPCSolver PROPERTIES PREFIX "")
install(TARGETS PyBoundednessMPCSolver DESTINATION "${PYTHON_SITELIB}/${MODULE_NAME}")

PYTHON_INSTALL_BUILD(${MODULE_NAME} __init__.py "${PYTHON_SITELIB}")

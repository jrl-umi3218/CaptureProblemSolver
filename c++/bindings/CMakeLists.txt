INCLUDE(../cmake/python.cmake)
FINDPYTHON()
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

set(MODULE_NAME "BalanceMPCSolver")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}/__init__.py COPYONLY)

set(PYBMS_SOURCES PyBalanceMPCSolver.cpp)

set(BMS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(PYBMS_HEADERS converters.h ${BMS_INCLUDE_DIR}/bms_api.h)

add_library(PyBalanceMPCSolver SHARED ${PYBMS_SOURCES} ${PYBMS_HEADERS})

target_link_libraries(PyBalanceMPCSolver BalanceMPCSolver boost_numpy)
TARGET_LINK_BOOST_PYTHON(PyBalanceMPCSolver)

set_target_properties(PyBalanceMPCSolver PROPERTIES PREFIX "")
install(TARGETS PyBalanceMPCSolver DESTINATION "${PYTHON_SITELIB}/${MODULE_NAME}")

PYTHON_INSTALL_BUILD(${MODULE_NAME} __init__.py "${PYTHON_SITELIB}")

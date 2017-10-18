set(BOOST_COMPONENTS unit_test_framework timer system)
search_for_boost()
add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${Boost_INCLUDE_DIRS})
if(WIN32)
  #This is one of the way to avoid link errors related to static variables in program_options
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif(WIN32)

macro(addUnitTest name)
  add_executable(${name} ${name}.cpp ${ARGN})
  target_link_libraries(${name} ${Boost_LIBRARIES} BalanceMPCSolver)
  add_test(${name}Unit ${name})
  # Adding a project configuration file (for MSVC only)
  GENERATE_MSVC_DOT_USER_FILE(${name})
endmacro(addUnitTest)

addUnitTest(MatrixComputationsTest)
addUnitTest(LinearConstraintsTest)

set(TESTSOURCES
    main.cpp
)

set(TESTHEADERS
)

add_executable(main ${TESTSOURCES} ${TESTHEADERS})
target_link_libraries(main PUBLIC BalanceMPCSolver)
GENERATE_MSVC_DOT_USER_FILE(main)


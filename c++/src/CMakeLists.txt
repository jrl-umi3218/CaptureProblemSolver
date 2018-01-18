set(CPS_SOURCES
  BoundenessConstraint.cpp
  CondensedOrthogonalMatrix.cpp
  Givens.cpp
  GivensSequence.cpp
  LeastSquare.cpp
  LinearConstraints.cpp
  Problem.cpp
  ProblemMatrices.cpp
  QRAlgorithms.cpp
  QuadraticObjective.cpp
  SQP.cpp)

set(CPS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(CPS_HEADERS
  ${CPS_INCLUDE_DIR}/cps/cps_api.h
  ${CPS_INCLUDE_DIR}/cps/BoundenessConstraint.h
  ${CPS_INCLUDE_DIR}/cps/CondensedOrthogonalMatrix.h
  ${CPS_INCLUDE_DIR}/cps/defs.h
  ${CPS_INCLUDE_DIR}/cps/Givens.h
  ${CPS_INCLUDE_DIR}/cps/GivensSequence.h
  ${CPS_INCLUDE_DIR}/cps/LeastSquare.h
  ${CPS_INCLUDE_DIR}/cps/LinearConstraints.h
  ${CPS_INCLUDE_DIR}/cps/Problem.h
  ${CPS_INCLUDE_DIR}/cps/ProblemMatrices.h
  ${CPS_INCLUDE_DIR}/cps/QRAlgorithms.h
  ${CPS_INCLUDE_DIR}/cps/QuadraticObjective.h
  ${CPS_INCLUDE_DIR}/cps/SQP.h
  ${CPS_INCLUDE_DIR}/cps/Statistics.h
  ${CPS_INCLUDE_DIR}/cps/toMatlab.h)
  
add_library(CaptureProblemSolver SHARED ${CPS_SOURCES} ${CPS_HEADERS})
set_target_properties(CaptureProblemSolver PROPERTIES COMPILE_FLAGS "-DCPS_EXPORTS")

if(STATISTICS)
  add_definitions(-DUSE_STATS)
endif()

install(TARGETS CaptureProblemSolver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CPS_INCLUDE_DIR}/ DESTINATION include
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

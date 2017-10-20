set(BMS_SOURCES
  BoundenessConstraint.cpp
  Givens.cpp
  GivensSequence.cpp
  LeastSquare.cpp
  LinearConstraints.cpp
  Logger.cpp
  Problem.cpp
  ProblemMatrices.cpp
  QRAlgorithms.cpp
  SQP.cpp)

set(BMS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(BMS_HEADERS
  ${BMS_INCLUDE_DIR}/bms_api.h
  ${BMS_INCLUDE_DIR}/BoundenessConstraint.h
  ${BMS_INCLUDE_DIR}/defs.h
  ${BMS_INCLUDE_DIR}/Givens.h
  ${BMS_INCLUDE_DIR}/GivensSequence.h
  ${BMS_INCLUDE_DIR}/LeastSquare.h
  ${BMS_INCLUDE_DIR}/LinearConstraints.h
  ${BMS_INCLUDE_DIR}/Logger.h
  ${BMS_INCLUDE_DIR}/Problem.h
  ${BMS_INCLUDE_DIR}/ProblemMatrices.h
  ${BMS_INCLUDE_DIR}/QRAlgorithms.h
  ${BMS_INCLUDE_DIR}/QuadraticObjective.h
  ${BMS_INCLUDE_DIR}/SQP.h
  ${BMS_INCLUDE_DIR}/toMatlab.h)
  
add_library(BalanceMPCSolver SHARED ${BMS_SOURCES} ${BMS_HEADERS})
set_target_properties(BalanceMPCSolver PROPERTIES COMPILE_FLAGS "-DBMS_EXPORTS")

install(TARGETS BalanceMPCSolver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${TVM_INCLUDE_DIR} DESTINATION include
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
set(BMS_SOURCES
  BoundenessConstraint.cpp
  CondensedOrthogonalMatrix.cpp
  Givens.cpp
  GivensSequence.cpp
  LeastSquare.cpp
  LinearConstraints.cpp
  Logger.cpp
  Problem.cpp
  ProblemMatrices.cpp
  QRAlgorithms.cpp
  QuadraticObjective.cpp
  SQP.cpp)

set(BMS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../include)
set(BMS_HEADERS
  ${BMS_INCLUDE_DIR}/bms/bms_api.h
  ${BMS_INCLUDE_DIR}/bms/BoundenessConstraint.h
  ${BMS_INCLUDE_DIR}/bms/CondensedOrthogonalMatrix.h
  ${BMS_INCLUDE_DIR}/bms/defs.h
  ${BMS_INCLUDE_DIR}/bms/Givens.h
  ${BMS_INCLUDE_DIR}/bms/GivensSequence.h
  ${BMS_INCLUDE_DIR}/bms/LeastSquare.h
  ${BMS_INCLUDE_DIR}/bms/LinearConstraints.h
  ${BMS_INCLUDE_DIR}/bms/Logger.h
  ${BMS_INCLUDE_DIR}/bms/Problem.h
  ${BMS_INCLUDE_DIR}/bms/ProblemMatrices.h
  ${BMS_INCLUDE_DIR}/bms/QRAlgorithms.h
  ${BMS_INCLUDE_DIR}/bms/QuadraticObjective.h
  ${BMS_INCLUDE_DIR}/bms/SQP.h
  ${BMS_INCLUDE_DIR}/bms/Statistics.h
  ${BMS_INCLUDE_DIR}/bms/toMatlab.h)
  
add_library(BoundednessMPCSolver SHARED ${BMS_SOURCES} ${BMS_HEADERS})
set_target_properties(BoundednessMPCSolver PROPERTIES COMPILE_FLAGS "-DBMS_EXPORTS")

if(STATISTICS)
  add_definitions(-DUSE_STATS)
endif()

install(TARGETS BoundednessMPCSolver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${BMS_INCLUDE_DIR}/ DESTINATION include
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

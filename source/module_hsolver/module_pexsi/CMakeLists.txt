add_library(pexsi OBJECT DistBCDMatrix.cpp DistCCSMatrix.cpp DistMatrixTransformer.cpp pexsi_solver.cpp simplePEXSI.cpp)

if(ENABLE_COVERAGE)
  add_coverage(pexsi)
endif()

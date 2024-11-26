# add diago shared library
list(APPEND _diago
    ${HSOLVER_PATH}/diago_dav_subspace.cpp
    ${HSOLVER_PATH}/diago_david.cpp
    ${HSOLVER_PATH}/diago_cg.cpp
    ${HSOLVER_PATH}/diag_const_nums.cpp
    ${HSOLVER_PATH}/diago_iter_assist.cpp

    ${HSOLVER_PATH}/kernels/dngvd_op.cpp
    ${HSOLVER_PATH}/kernels/math_kernel_op.cpp
    # dependency
    ${BASE_PATH}/kernels/math_op.cpp
    ${BASE_PATH}/module_device/device.cpp
    ${BASE_PATH}/module_device/memory_op.cpp
    
    ${PSI_PATH}/psi.cpp
  )
add_library(diagopack SHARED 
    ${_diago}
    )
target_link_libraries(diagopack
    base
    parameter
    container
    orb
    ${math_libs}
    ${OpenBLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
    )

list(APPEND pymodule_hsolver
    ${PROJECT_SOURCE_DIR}/src/hsolver/py_hsolver.cpp
    )

pybind11_add_module(_hsolver_pack MODULE ${pymodule_hsolver})
target_link_libraries(_hsolver_pack PRIVATE pybind11::headers diagopack)
target_compile_definitions(_hsolver_pack PRIVATE VERSION_INFO=${PROJECT_VERSION})

set_target_properties(diagopack PROPERTIES INSTALL_RPATH "$ORIGIN")
set_target_properties(_hsolver_pack PROPERTIES INSTALL_RPATH "$ORIGIN")

install(TARGETS _hsolver_pack diagopack DESTINATION ${TARGET_PACK}/hsolver)
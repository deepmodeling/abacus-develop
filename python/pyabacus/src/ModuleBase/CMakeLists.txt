list(APPEND pymodule_base
    ${PROJECT_SOURCE_DIR}/src/ModuleBase/py_base_math.cpp
    ${BASE_PATH}/kernels/math_op.cpp
    ${BASE_PATH}/kernels/math_kernel_op.cpp
    ${BASE_PATH}/module_device/memory_op.cpp
    ${BASE_PATH}/module_device/device.cpp
    )

pybind11_add_module(_base_pack MODULE ${pymodule_base})

target_link_libraries(_base_pack
    PRIVATE
    base
    parameter
    container
    orb
    ${math_libs}
    ${OpenBLAS_LIBRARIES} 
    ${LAPACK_LIBRARIES}
    )

target_link_libraries(_base_pack PRIVATE pybind11::headers)
target_compile_definitions(_base_pack PRIVATE VERSION_INFO=${PROJECT_VERSION})

set_target_properties(_base_pack PROPERTIES INSTALL_RPATH "$ORIGIN")

install(TARGETS _base_pack DESTINATION ${TARGET_PACK}/ModuleBase)
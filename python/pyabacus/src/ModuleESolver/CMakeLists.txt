# CMakeLists.txt for ModuleESolver Python bindings
# This module provides Python bindings for ESolver_KS_LCAO

# Set paths for ESolver related sources
set(ESOLVER_PATH "${ABACUS_SOURCE_DIR}/source_esolver")
set(ESTATE_PATH "${ABACUS_SOURCE_DIR}/source_estate")
set(LCAO_PATH "${ABACUS_SOURCE_DIR}/source_lcao")
set(HAMILT_PATH "${ABACUS_SOURCE_DIR}/source_hamilt")
set(CELL_PATH "${ABACUS_SOURCE_DIR}/source_cell")
set(PSI_PATH "${ABACUS_SOURCE_DIR}/source_psi")
set(BASIS_PATH "${ABACUS_SOURCE_DIR}/source_basis")
set(IO_PATH "${ABACUS_SOURCE_DIR}/source_io")

# Python module source files - only the binding code
list(APPEND pymodule_esolver
    ${PROJECT_SOURCE_DIR}/src/ModuleESolver/py_esolver_lcao.cpp
)

# Create pybind11 module
pybind11_add_module(_esolver_pack MODULE ${pymodule_esolver})

target_include_directories(_esolver_pack PRIVATE
    ${ABACUS_SOURCE_DIR}
    ${ESOLVER_PATH}
    ${ESTATE_PATH}
    ${ESTATE_PATH}/module_charge
    ${ESTATE_PATH}/module_dm
    ${ESTATE_PATH}/potentials
    ${LCAO_PATH}
    ${LCAO_PATH}/module_hcontainer
    ${HAMILT_PATH}
    ${CELL_PATH}
    ${PSI_PATH}
    ${BASIS_PATH}
    ${BASIS_PATH}/module_ao
    ${BASIS_PATH}/module_nao
    ${IO_PATH}
    ${IO_PATH}/module_parameter
)

# Only link pybind11 headers - the module uses placeholder implementations
# that don't require actual ABACUS libraries
target_link_libraries(_esolver_pack PRIVATE
    pybind11::headers
)

target_compile_definitions(_esolver_pack PRIVATE VERSION_INFO=${PROJECT_VERSION})

# Set RPATH for shared libraries
set_target_properties(_esolver_pack PROPERTIES INSTALL_RPATH "$ORIGIN")

# Install targets
install(TARGETS _esolver_pack DESTINATION ${TARGET_PACK}/esolver)

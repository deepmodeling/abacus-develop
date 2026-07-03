file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/support DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(DEEPKS_UT_PP_ORB_DIR "${PROJECT_SOURCE_DIR}/tests/PP_ORB")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/support/NO_GO_deepks_UT/STRU
    ${CMAKE_CURRENT_BINARY_DIR}/support/NO_GO_deepks_UT/STRU
    @ONLY
)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/support/NO_KP_deepks_UT/STRU
    ${CMAKE_CURRENT_BINARY_DIR}/support/NO_KP_deepks_UT/STRU
    @ONLY
)

set(DEEPKS_UNIT_COMMON_SOURCES
    deepks_test_prep.cpp
    deepks_test.cpp
    ../../../source_cell/unitcell.cpp
    ../../../source_cell/update_cell.cpp
    ../../../source_cell/bcast_cell.cpp
    ../../../source_cell/atom_spec.cpp
    ../../../source_cell/atom_pseudo.cpp
    ../../../source_cell/read_atoms.cpp
    ../../../source_cell/read_atoms_helper.cpp
    ../../../source_cell/read_stru.cpp
    ../../../source_cell/print_cell.cpp
    ../../../source_cell/read_atom_species.cpp
    ../../../source_cell/klist.cpp
    ../../../source_cell/parallel_kpoints.cpp
    ../../../source_cell/k_vector_utils.cpp
    ../../../source_cell/setup_nonlocal.cpp
    ../../../source_cell/pseudo.cpp
    ../../../source_cell/read_pp.cpp
    ../../../source_cell/read_pp_complete.cpp
    ../../../source_cell/read_pp_upf100.cpp
    ../../../source_cell/read_pp_upf201.cpp
    ../../../source_cell/read_pp_vwr.cpp
    ../../../source_cell/read_pp_blps.cpp
    ../../../source_cell/sep.cpp
    ../../../source_cell/sep_cell.cpp
    ../../../source_pw/module_pwdft/soc.cpp
    ../../../source_io/module_output/sparse_matrix.cpp
    ../../../source_estate/read_pseudo.cpp
    ../../../source_estate/cal_wfc.cpp
    ../../../source_estate/read_orb.cpp
    ../../../source_estate/cal_nelec_nband.cpp
    ../../../source_estate/module_dm/density_matrix.cpp
    ../../../source_estate/module_dm/density_matrix_io.cpp
    ../../center2_orb.cpp
    ../../center2_orb-orb11.cpp
    ../../center2_orb-orb21.cpp
    ../../module_hcontainer/base_matrix.cpp
    ../../module_hcontainer/hcontainer.cpp
    ../../module_hcontainer/atom_pair.cpp
    ../../module_hcontainer/func_transfer.cpp
    ../../module_hcontainer/func_folding.cpp
    ../../module_hcontainer/transfer.cpp
    ../../module_hcontainer/output_hcontainer.cpp
    ../../module_operator_lcao/deepks_lcao.cpp
    ../../module_operator_lcao/operator_lcao.cpp
    ../../../source_hamilt/operator.cpp
    ../../../source_io/module_hs/cal_r_overlap_R.cpp
    ../../../source_io/module_hs/single_R_io.cpp
    ../../../source_io/module_hs/rr_sparse_writer.cpp
    ../../module_rt/td_folding.cpp
    mock_berryphase.cpp
    mock_tdinfo.cpp
)

add_library(deepks_unit_support OBJECT ${DEEPKS_UNIT_COMMON_SOURCES})

if(ENABLE_COVERAGE)
  add_coverage(deepks_unit_support)
endif()

set(DEEPKS_UNIT_LIBS
    deepks_unit_support
    base
    device
    parameter
    deepks
    psi
    planewave
    neighbor
    container
    orb
    gint
    numerical_atomic_orbitals
    symmetry
    ${math_libs}
)

set(DEEPKS_UNIT_PHIALPHA_SOURCES
    deepks_test_phialpha.cpp
)
set(DEEPKS_UNIT_PDM_SOURCES
    ${DEEPKS_UNIT_PHIALPHA_SOURCES}
    deepks_test_pdm.cpp
)
set(DEEPKS_UNIT_DESCRIPTOR_SOURCES
    ${DEEPKS_UNIT_PDM_SOURCES}
    deepks_test_descriptor.cpp
)
set(DEEPKS_UNIT_GDMX_SOURCES
    ${DEEPKS_UNIT_PDM_SOURCES}
    deepks_test_gdmx.cpp
)
set(DEEPKS_UNIT_GVX_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_gdmx.cpp
    deepks_test_gvx.cpp
)
set(DEEPKS_UNIT_GDMEPSL_SOURCES
    ${DEEPKS_UNIT_PDM_SOURCES}
    deepks_test_gdmepsl.cpp
)
set(DEEPKS_UNIT_GVEPSL_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_gdmepsl.cpp
    deepks_test_gvepsl.cpp
)
set(DEEPKS_UNIT_ORBPRE_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_orbpre.cpp
)
set(DEEPKS_UNIT_VDPRE_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_vdpre.cpp
)
set(DEEPKS_UNIT_VDRPRE_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_vdrpre.cpp
)
set(DEEPKS_UNIT_EDELTA_SOURCES
    ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
    deepks_test_edelta.cpp
)
set(DEEPKS_UNIT_E_DELTABANDS_SOURCES
    ${DEEPKS_UNIT_EDELTA_SOURCES}
    deepks_test_e_deltabands.cpp
)
set(DEEPKS_UNIT_FORCE_STRESS_DELTA_SOURCES
    ${DEEPKS_UNIT_EDELTA_SOURCES}
    deepks_test_force_stress_delta.cpp
)
set(DEEPKS_UNIT_O_DELTA_SOURCES
    ${DEEPKS_UNIT_E_DELTABANDS_SOURCES}
    deepks_test_o_delta.cpp
)

function(configure_deepks_unit_target TARGET_NAME CHECK_NAME CASE_DIR)
  target_compile_definitions(
      ${TARGET_NAME}
      PRIVATE
      DEEPKS_UT_CHECK_NAME="${CHECK_NAME}"
      DEEPKS_UT_CASE_DIR="${CASE_DIR}"
      DEEPKS_UT_RUNNER=run_deepks_unit_${CHECK_NAME}
  )
endfunction()

AddTest(
    TARGET MODULE_LCAO_DEEPKS_phialpha_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_PHIALPHA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_phialpha_gamma phialpha NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_pdm_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_PDM_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_pdm_gamma pdm NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_descriptor_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_descriptor_gamma descriptor NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gdmx_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GDMX_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gdmx_gamma gdmx NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gvx_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GVX_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gvx_gamma gvx NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gdmepsl_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GDMEPSL_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gdmepsl_gamma gdmepsl NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gvepsl_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GVEPSL_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gvepsl_gamma gvepsl NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_orbpre_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_ORBPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_orbpre_gamma orbpre NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_vdpre_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_VDPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_vdpre_gamma vdpre NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_vdrpre_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_VDRPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_vdrpre_gamma vdrpre NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_edelta_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_EDELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_edelta_gamma edelta NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_e_deltabands_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_E_DELTABANDS_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_e_deltabands_gamma e_deltabands NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_force_stress_delta_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_FORCE_STRESS_DELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_force_stress_delta_gamma force_stress_delta NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_o_delta_gamma
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_O_DELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_o_delta_gamma o_delta NO_GO_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_phialpha_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_PHIALPHA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_phialpha_multik phialpha NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_pdm_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_PDM_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_pdm_multik pdm NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_descriptor_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_DESCRIPTOR_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_descriptor_multik descriptor NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gdmx_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GDMX_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gdmx_multik gdmx NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gvx_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GVX_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gvx_multik gvx NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gdmepsl_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GDMEPSL_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gdmepsl_multik gdmepsl NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_gvepsl_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_GVEPSL_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_gvepsl_multik gvepsl NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_orbpre_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_ORBPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_orbpre_multik orbpre NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_vdpre_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_VDPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_vdpre_multik vdpre NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_vdrpre_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_VDRPRE_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_vdrpre_multik vdrpre NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_edelta_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_EDELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_edelta_multik edelta NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_e_deltabands_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_E_DELTABANDS_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_e_deltabands_multik e_deltabands NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_force_stress_delta_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_FORCE_STRESS_DELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_force_stress_delta_multik force_stress_delta NO_KP_deepks_UT)

AddTest(
    TARGET MODULE_LCAO_DEEPKS_o_delta_multik
    LIBS ${DEEPKS_UNIT_LIBS}
    SOURCES main_deepks.cpp ${DEEPKS_UNIT_O_DELTA_SOURCES}
)
configure_deepks_unit_target(MODULE_LCAO_DEEPKS_o_delta_multik o_delta NO_KP_deepks_UT)

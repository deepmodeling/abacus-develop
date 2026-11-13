abacus_disable_feature_definitions(__CUDA)

if(ENABLE_LCAO)

  addtest(
    TARGET MODULE_LCAO_deltaspin_template_helpers
    LIBS base device parameter symmetry
    SOURCES test_template_helpers.cpp
            ../spin_constrain.cpp
            ../deltaspin_state.cpp
            ../template_helpers.cpp
            ../lambda_loop_helper.cpp
            ../basic_funcs.cpp
            ../../../source_cell/klist.cpp
            ../../../source_cell/klist_io.cpp
            ../../../source_cell/parallel_kpoints.cpp
            ../../../source_cell/reciprocal_grid.cpp)

endif()

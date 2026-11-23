abacus_disable_feature_definitions(__MPI __LCAO)

addtest(
  TARGET MODULE_ESOLVER_forcefields
  LIBS parameter base device
  SOURCES test_esolver_forcefields.cpp ../esolver_forcefields.cpp
          ../../source_cell/basecell.cpp
          ../../source_cell/atom_spec.cpp
          ../../source_cell/atom_pseudo.cpp
          ../../source_cell/pseudo.cpp
          ../../source_cell/unitcell.cpp
          ../../source_cell/magnetism.cpp
          ../../source_cell/mdcell.cpp
          ../../source_cell/sep_cell.cpp
          ../../source_cell/sep.cpp
          ../../source_io/module_output/output_log.cpp
          ../../source_io/module_output/print_info.cpp)

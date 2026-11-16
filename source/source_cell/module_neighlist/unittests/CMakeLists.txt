abacus_disable_feature_definitions(__MPI)
abacus_add_local_feature_definitions(__NORMAL)

addtest(
  TARGET MODULE_CELL_NEIGHBOR_domain_decomposition_serial
  LIBS parameter base device
  SOURCES test_domain_decomposition.cpp
          ../../mdcell.cpp
          ../../dom_core.cpp
          ../../dom_ghost.cpp
          ../../dom_migrate.cpp
          ../../dom_force.cpp
          ../neighbor_search.cpp
          ../bin_manager.cpp
          ../page_allocator.cpp)

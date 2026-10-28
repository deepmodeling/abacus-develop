if(ENABLE_LCAO)
  add_library(
      orb
      OBJECT
      orb_atomic.cpp
      orb_atomic_lm.cpp
      orb_gaunt_table.cpp
      orb_nonlocal.cpp
      orb_nonlocal_lm.cpp
      orb_read.cpp
      parallel_orbitals.cpp
      elem_basis_idx_orb.cpp
  )

  if(ENABLE_COVERAGE)
    add_coverage(orb)
  endif()

  if(BUILD_TESTING)
    # module_ao 测试主体均以 #ifdef __MPI 守卫或纯串行，支持 MPI/串行双配置
    # （D6：set_serial 与 MPI 双覆盖，串行亦须构建并运行本目录测试）
    add_subdirectory(test)
  endif()
endif()

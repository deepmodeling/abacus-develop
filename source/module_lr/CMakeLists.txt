if(ENABLE_LCAO)
  add_subdirectory(utils)
  add_subdirectory(ao_to_mo_transformer)
  add_subdirectory(dm_trans)
  add_subdirectory(ri_benchmark)

  list(APPEND objects
    utils/lr_util.cpp
    utils/lr_util_hcontainer.cpp
    ao_to_mo_transformer/ao_to_mo_parallel.cpp
    ao_to_mo_transformer/ao_to_mo_serial.cpp
    dm_trans/dm_trans_parallel.cpp
    dm_trans/dm_trans_serial.cpp
    dm_trans/dmr_complex.cpp
    operator_casida/operator_lr_hxc.cpp
    operator_casida/operator_lr_exx.cpp
    potentials/pot_hxc_lrtd.cpp
    lr_spectrum.cpp
    lr_spectrum_velocity.cpp
    hamilt_casida.cpp
    esolver_lrtd_lcao.cpp
    potentials/xc_kernel.cpp)

add_library(
    lr
    OBJECT
    ${objects}
)

endif()
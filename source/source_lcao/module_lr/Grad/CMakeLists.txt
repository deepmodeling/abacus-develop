add_subdirectory(dm_diff)
add_subdirectory(CVCX)

add_library(
lr_grad
OBJECT
CVCX/CVCX_parallel.cpp
CVCX/CVCX_serial.cpp
xc/pot_grad_xc.cpp
force/lr_force.cpp
force/lr_force_test.cpp
multipliers/cal_edm_from_multipliers.cpp
esolver_lr_grad.cpp
)
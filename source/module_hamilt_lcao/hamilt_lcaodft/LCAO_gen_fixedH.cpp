#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include <vector>
#include <unordered_map>
#include <map>
#include "module_base/timer.h"

#ifdef __MKL
#include <mkl_service.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

LCAO_gen_fixedH::LCAO_gen_fixedH()
{}

LCAO_gen_fixedH::~LCAO_gen_fixedH()
{}


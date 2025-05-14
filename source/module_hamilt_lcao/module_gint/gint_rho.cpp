#include "gint_k.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/array_pool.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint::cal_meshball_rho(const int na_grid,
                            const int*const block_index,
                            const int*const vindex,
                            const double*const*const psir_ylm,
                            const double*const*const psir_DMR,
                            double*const rho)
{
    const int inc = 1;
    // sum over mu to get density on grid
    for (int ib = 0; ib < this->bxyz; ++ib)
    {
        const double r = ddot_(&block_index[na_grid], psir_ylm[ib], &inc, psir_DMR[ib], &inc);
        const int grid = vindex[ib];
        rho[grid] += r;
    }
}

#ifndef HSOLVER_H
#define HSOLVER_H

#include "module_base/macros.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_general/hamilt.h"
#include "module_psi/wavefunc.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_psi/psi.h"

#include <complex>

namespace hsolver
{

double set_diagethr_ks(const std::string basis_type,
                       const std::string esolver_type,
                       const std::string calculation_in,
                       const std::string init_chg_in,
                       const std::string precision_flag_in,
                       const int istep,
                       const int iter,
                       const double drho,
                       const double pw_diag_thr_init,
                       const double diag_ethr_in,
                       const double nelec_in);

double set_diagethr_sdft(const std::string basis_type,
                         const std::string esolver_type,
                         const std::string calculation_in,
                         const std::string init_chg_in,
                         const int istep,
                         const int iter,
                         const double drho,
                         const double pw_diag_thr_init,
                         const double diag_ethr_in,
                         const int nband_in,
                         const double stoiter_ks_ne_in);


// reset diagethr according to drho and hsolver_error
double reset_diag_ethr(std::ofstream& ofs_running,
                       const std::string basis_type,
                       const std::string esolver_type,
                       const std::string precision_flag_in,
                       const double hsover_error,
                       const double drho_in,
                       const double diag_ethr_in,
                       const double nelec_in);

// calculate hsolver_error (for sdft, lcao and lcao-in-pw, we suppose the error is zero)
double cal_hsolve_error(const std::string basis_type,
                        const std::string esolver_type,
                        const double diag_ethr_in,
                        const double nelec_in);

} // namespace hsolver
#endif
#ifndef HSOLVER_H
#define HSOLVER_H

#include "diagh.h"
#include "module_base/macros.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_psi/psi.h"

#include <complex>

namespace hsolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolver
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    HSolver() {};

    // solve Hamiltonian to electronic density in ElecState
    virtual void solve(hamilt::Hamilt<T, Device>* phm,
                       psi::Psi<T, Device>& ppsi,
                       elecstate::ElecState* pes,
                       const std::string method,
                       const bool skip_charge)
    {
        return;
    }

    virtual void solve(hamilt::Hamilt<T, Device>* phm,
                       psi::Psi<T, Device>& ppsi,
                       elecstate::ElecState* pes,
                       ModulePW::PW_Basis_K* wfc_basis,
                       Stochastic_WF& stowf,
                       const int istep,
                       const int iter,
                       const std::string method,
                       const int scf_iter_in,
                       const bool need_subspace_in,
                       const int diag_iter_max_in,
                       const double pw_diag_thr_in,
                       const bool skip_charge)
    {
        return;
    }

    // set diagethr according to drho (for lcao and lcao-in-pw, we suppose the error is zero and we set diagethr to 0)
    virtual Real set_diagethr(Real diag_ethr_in, const int istep, const int iter, const Real drho)
    {
        return 0.0;
    }

    // reset diagethr according to drho and hsolver_error
    virtual Real reset_diagethr(std::ofstream& ofs_running, const Real hsover_error, const Real drho, Real diag_ethr_in)
    {
        return 0.0;
    }

    // calculate hsolver_error (for sdft, lcao and lcao-in-pw, we suppose the error is zero)
    virtual Real cal_hsolerror(const Real diag_ethr_in)
    {
        return 0.0;
    };
};

double reset_diag_ethr(std::ofstream& ofs_running,
                       const double hsover_error,
                       const double drho,
                       double diag_ethr_in,
                       std::string basis_type,
                       std::string esolver_type);

double cal_hsolve_error(const double diag_ethr_in, std::string basis_type, std::string esolver_type);

} // namespace hsolver
#endif
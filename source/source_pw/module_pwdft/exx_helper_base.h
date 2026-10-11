#ifndef EXX_HELPER_BASE_H
#define EXX_HELPER_BASE_H

#include "source_base/matrix.h"

class Charge;
class UnitCell;
struct Input_para;
struct General_Exx_Info;

class Exx_HelperBase
{
  public:
    Exx_HelperBase() = default;
    virtual ~Exx_HelperBase() = default;

    virtual void init(const UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& wg, const General_Exx_Info& exx_info) = 0;

    virtual void before_scf(void* p_hamilt, void* psi, const Input_para& inp, const General_Exx_Info& exx_info) = 0;

    virtual bool iter_finish(void* p_elec, Charge* p_charge, void* psi,
                             UnitCell& ucell, const Input_para& inp,
                             bool& conv_esolver, int iter) = 0;

    virtual void set_firstiter(bool flag = true) = 0;
    virtual void set_wg(const ModuleBase::matrix* wg) = 0;
    virtual void set_psi(void* psi) = 0;
    virtual void iter_inc() = 0;

    virtual void set_op() = 0;

    /// @brief Decide whether EXX has converged and SCF should rerun.
    /// @note iter is taken by value: EXX no longer mutates the caller's
    ///       iteration counter. SCF rerun is requested via the return
    ///       value (false = not converged, rerun SCF) and the caller
    ///       (ESolver_KS_PW) sets ESolver_KS::scf_rerun_ to trigger the
    ///       restart at iter=1 in ESolver_KS::runner.
    virtual bool exx_after_converge(int iter, bool ene_conv) = 0;

    virtual double cal_exx_energy(void* psi) = 0;

    virtual bool get_op_first_iter() const = 0;
    virtual void set_op_first_iter(bool flag) = 0;

    virtual void set_op_exx(void* op) = 0;
};

#endif // EXX_HELPER_BASE_H

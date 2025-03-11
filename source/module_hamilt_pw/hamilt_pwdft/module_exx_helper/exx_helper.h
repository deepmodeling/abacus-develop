//
// For EXX in PW.
//
#include "module_psi/psi.h"
#include "module_base/matrix.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/op_exx_pw.h"

#ifndef EXX_HELPER_H
#define EXX_HELPER_H
template <typename T, typename Device>
struct Exx_Helper
{
    using Real = typename GetTypeReal<T>::type;
    using OperatorEXX = hamilt::OperatorEXXPW<T, Device>;

  public:
    Exx_Helper() = default;
    OperatorEXX *op_exx = nullptr;

    void set_firstiter() { first_iter = true; }
    void set_wg(const ModuleBase::matrix *wg_) { wg = wg_; }
    void set_psi(psi::Psi<T, Device> *psi_);

    void set_op()
    {
        op_exx->first_iter = first_iter;
        set_psi(psi);
        op_exx->set_wg(wg);
    }

    bool exx_after_converge(int &iter);

    double cal_exx_energy(psi::Psi<T, Device> *psi_);

  private:
    bool first_iter;
    psi::Psi<T, Device> *psi = nullptr;
    const ModuleBase::matrix *wg = nullptr;

};
#endif // EXX_HELPER_H

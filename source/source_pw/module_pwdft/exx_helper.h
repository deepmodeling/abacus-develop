#include "source_psi/psi.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/op_pw_exx.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_pw/module_pwdft/exx_helper_base.h"

#ifndef EXX_HELPER_H
#define EXX_HELPER_H

class Charge;

template <typename T, typename Device>
struct Exx_Helper : public Exx_HelperBase
{
    using Real = typename GetTypeReal<T>::type;
    using OperatorEXX = hamilt::OperatorEXXPW<T, Device>;

  public:
    Exx_Helper() = default;
    virtual ~Exx_Helper() = default;
    OperatorEXX *op_exx = nullptr;

    void init(const UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& wg) override;

    void before_scf(void* p_hamilt, void* psi, const Input_para& inp) override;

    bool iter_finish(void* p_elec, Charge* p_charge, void* psi,
                     UnitCell& ucell, const Input_para& inp,
                     bool& conv_esolver, int& iter) override;

    void set_firstiter(bool flag = true) override { first_iter = flag; }
    void set_wg(const ModuleBase::matrix *wg_) override { wg = wg_; }
    void set_psi(void* psi_) override;
    void iter_inc() override { exx_iter++; }

    void set_op() override
    {
        op_exx->first_iter = first_iter;
        set_psi(psi);
        op_exx->set_wg(wg);
    }

    bool exx_after_converge(int &iter, bool ene_conv) override;

    double cal_exx_energy(void* psi_) override;

    bool get_op_first_iter() const override { return op_exx ? op_exx->first_iter : false; }
    void set_op_first_iter(bool flag) override { if (op_exx) op_exx->first_iter = flag; }
    void set_op_exx(void* op) override { op_exx = reinterpret_cast<OperatorEXX*>(op); }

  private:
    bool first_iter = false;
    psi::Psi<T, Device> *psi = nullptr;
    const ModuleBase::matrix *wg = nullptr;
    int exx_iter = 0;

};
#endif // EXX_HELPER_H

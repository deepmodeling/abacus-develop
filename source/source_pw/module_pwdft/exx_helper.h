#include "source_psi/psi.h"
#include "source_base/matrix.h"
#include "source_pw/module_pwdft/op_pw_exx.h"
#include "source_io/module_parameter/input_parameter.h"

#ifndef EXX_HELPER_H
#define EXX_HELPER_H

class Charge;

template <typename T, typename Device>
struct Exx_Helper
{
    using Real = typename GetTypeReal<T>::type;
    using OperatorEXX = hamilt::OperatorEXXPW<T, Device>;

  public:
    Exx_Helper() = default;
    OperatorEXX *op_exx = nullptr;

    void init(const UnitCell& ucell, const Input_para& inp, const ModuleBase::matrix& wg);

    /**
     * @brief Setup EXX helper before SCF iteration.
     *
     * This function sets up the EXX helper for the Hamiltonian and psi
     * before each SCF iteration. It checks if the calculation type and
     * EXX settings are appropriate.
     *
     * @param p_hamilt Pointer to the Hamiltonian object (void* to avoid circular dependency).
     * @param psi Pointer to the wave function object.
     * @param inp The input parameters.
     */
    void before_scf(void* p_hamilt, psi::Psi<T, Device>* psi, const Input_para& inp);

    /**
     * @brief Handle EXX-related operations after SCF iteration.
     *
     * This function handles EXX convergence checking and potential update
     * after each SCF iteration. It is called in iter_finish.
     *
     * @param p_elec Pointer to the ElecState object (void* to avoid circular dependency).
     * @param p_charge Pointer to the Charge object.
     * @param psi Pointer to the wave function object.
     * @param ucell The unit cell (non-const reference for update_pot).
     * @param inp The input parameters.
     * @param conv_esolver Whether SCF is converged (may be modified).
     * @param iter The current iteration number (may be modified).
     * @return true if EXX processing was done, false otherwise.
     */
    bool iter_finish(void* p_elec, Charge* p_charge, psi::Psi<T, Device>* psi,
                     UnitCell& ucell, const Input_para& inp,
                     bool& conv_esolver, int& iter);

    void set_firstiter(bool flag = true) { first_iter = flag; }
    void set_wg(const ModuleBase::matrix *wg_) { wg = wg_; }
    void set_psi(psi::Psi<T, Device> *psi_);
    void iter_inc() { exx_iter++; }

    void set_op()
    {
        op_exx->first_iter = first_iter;
        set_psi(psi);
        op_exx->set_wg(wg);
    }

    bool exx_after_converge(int &iter, bool ene_conv);

    double cal_exx_energy(psi::Psi<T, Device> *psi_);

  private:
    bool first_iter = false;
    psi::Psi<T, Device> *psi = nullptr;
    const ModuleBase::matrix *wg = nullptr;
    int exx_iter = 0;

};
#endif // EXX_HELPER_H

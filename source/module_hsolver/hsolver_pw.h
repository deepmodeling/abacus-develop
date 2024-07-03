#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_base/macros.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"

namespace hsolver {

template <typename T, typename Device = base_device::DEVICE_CPU>
class HSolverPW : public HSolver<T, Device> {
  private:
    bool is_first_scf = true;

    // Note GetTypeReal<T>::type will
    // return T if T is real type(float, double),
    // otherwise return the real type of T(complex<float>, complex<double>)
    using Real = typename GetTypeReal<T>::type;

  public:
    /**
     * @brief diago_full_acc
     * If .TRUE. all the empty states are diagonalized at the same level of
     * accuracy of the occupied ones. Otherwise the empty states are
     * diagonalized using a larger threshold (this should not affect total
     * energy, forces, and other ground-state properties).
     *
     */
    static bool diago_full_acc;

    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in, wavefunc* pwf_in);

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    /// @brief solve function for pw
    /// @param pHamilt interface to hamilt
    /// @param psi reference to psi
    /// @param pes interface to elecstate
    /// @param method_in dav or cg
    /// @param skip_charge
    void solve(hamilt::Hamilt<T, Device>* pHamilt,
               psi::Psi<T, Device>& psi,
               elecstate::ElecState* pes,
               const std::string method_in,
               const bool skip_charge) override;
    /// @brief solve function for lcao_in_pw
    /// @param pHamilt interface to hamilt
    /// @param psi reference to psi
    /// @param pes interface to elecstate
    /// @param transform transformation matrix between lcao and pw
    /// @param skip_charge
    void solve(hamilt::Hamilt<T, Device>* pHamilt,
               psi::Psi<T, Device>& psi,
               elecstate::ElecState* pes,
               psi::Psi<T, Device>& transform,
               const bool skip_charge) override;
    virtual Real cal_hsolerror() override;
    virtual Real
        set_diagethr(const int istep, const int iter, const Real drho) override;
    virtual Real reset_diagethr(std::ofstream& ofs_running,
                                const Real hsover_error,
                                const Real drho) override;

  protected:
    // void initDiagh(const psi::Psi<T, Device>& psi_in);
    void endDiagh();
    void hamiltSolvePsiK(hamilt::Hamilt<T, Device>* hm,
                         psi::Psi<T, Device>& psi,
                         Real* eigenvalue);

    void updatePsiK(hamilt::Hamilt<T, Device>* pHamilt,
                    psi::Psi<T, Device>& psi,
                    const int ik);

    ModulePW::PW_Basis_K* wfc_basis = nullptr;
    wavefunc* pwf = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<Real>& h_diag,
                             const int ik,
                             const int npw);

    std::vector<Real> precondition;
    std::vector<Real> eigenvalues;
    std::vector<bool> is_occupied;

    bool initialed_psi = false;

    hamilt::Hamilt<T, Device>* hamilt_ = nullptr;

    Device* ctx = {};
    using resmem_var_op
        = base_device::memory::resize_memory_op<Real, base_device::DEVICE_CPU>;
    using delmem_var_op
        = base_device::memory::delete_memory_op<Real, base_device::DEVICE_CPU>;
    using castmem_2d_2h_op
        = base_device::memory::cast_memory_op<double,
                                              Real,
                                              base_device::DEVICE_CPU,
                                              base_device::DEVICE_CPU>;
};

template <typename T, typename Device>
bool HSolverPW<T, Device>::diago_full_acc = false;

} // namespace hsolver

#endif
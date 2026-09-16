#ifndef HAMILT_HS_ADAPTER_H
#define HAMILT_HS_ADAPTER_H

#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_hamilt/hamilt.h"
#include "source_hsolver/hs_matrix.h"
#include "source_hsolver/hs_operator.h"
#include "source_psi/psi.h"

namespace hamilt
{

/**
 * @brief Presents a Hamilt as the H/S block-vector operator the iterative
 *        eigensolvers work on (hsolver::HSOperator).
 *
 * This is the only place that wraps raw pointers into Psi / hpsi_info for the
 * operator chain; hsolver itself never sees Hamilt.
 */
template <typename T, typename Device = base_device::DEVICE_CPU>
class HamiltHSOperator : public hsolver::HSOperator<T, Device>
{
  public:
    HamiltHSOperator(Hamilt<T, Device>* hm, const ModulePW::PW_Basis_K* wfc_basis) : hm_(hm), wfc_basis_(wfc_basis)
    {
    }

    void update_k(const int ik) override
    {
        hm_->updateHk(ik);
        ik_ = ik;
        npw_ = wfc_basis_->npwk[ik];
    }

    void hpsi(const T* x, T* hx, const int ld, const int nvec) const override
    {
        if (hm_->ops == nullptr)
        {
            ModuleBase::WARNING_QUIT("HamiltHSOperator::hpsi", "Operators in Hamilt are not allocated yet");
        }
        // non-owning view of x: one k point, nvec bands, leading dimension ld, npw valid rows
        psi::Psi<T, Device> x_view(const_cast<T*>(x), 1, nvec, ld, npw_);
        typename Operator<T, Device>::hpsi_info info(&x_view, psi::Range(true, 0, 0, nvec - 1), hx);
        hm_->ops->hPsi(info);
    }

    void spsi(const T* x, T* sx, const int ld, const int nvec) const override
    {
        hm_->sPsi(x, sx, ld, npw_, nvec);
    }

  protected:
    Hamilt<T, Device>* hm_ = nullptr;
    const ModulePW::PW_Basis_K* wfc_basis_ = nullptr;
    int ik_ = 0;  ///< k point set by the last update_k()
    int npw_ = 0; ///< number of plane waves of that k point (without npol)
};

/**
 * @brief Presents a Hamilt as the H(k)/S(k) matrix source the direct
 *        eigensolvers work on (hsolver::HSMatrix).
 */
template <typename T>
class HamiltHSMatrix : public hsolver::HSMatrix<T>
{
  public:
    explicit HamiltHSMatrix(Hamilt<T>* hm) : hm_(hm)
    {
    }

    void hs_at_k(const int ik, ModuleBase::MatrixBlock<T>& hk, ModuleBase::MatrixBlock<T>& sk) override
    {
        hm_->updateHk(ik);
        hm_->matrix(hk, sk);
    }

  private:
    Hamilt<T>* hm_ = nullptr;
};

} // namespace hamilt

#endif // HAMILT_HS_ADAPTER_H

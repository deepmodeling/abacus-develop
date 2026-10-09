#ifndef PW_NONLOCAL_WORKSPACE_H
#define PW_NONLOCAL_WORKSPACE_H

#include "source_pw/module_pwdft/vnl_pw.h"
#include <cstddef>

namespace hamilt
{
/** @brief Reusable device buffers for applying a nonlocal projector operator. */
template <typename T, typename Device>
class NonlocalWorkspace
{
  private:
    T* becp_ = nullptr;
    T* ps_ = nullptr;
    size_t becp_capacity_ = 0;
    size_t ps_capacity_ = 0;

    void project(const T* vkb, const T* psi, const int npw, const int ldv,
                 const int ldp, const int nkb, const int bands);
    void contract(const UnitCell& cell, const pseudopot_cell_vnl& pp,
                  const int spin, const int npol, const int bands);
    void back_project(const T* vkb, T* hpsi, const int npw, const int ldv,
                      const int ldp, const int nkb, const int bands);

  public:
    NonlocalWorkspace() = default;
    ~NonlocalWorkspace();
    NonlocalWorkspace(const NonlocalWorkspace&) = delete;
    NonlocalWorkspace& operator=(const NonlocalWorkspace&) = delete;

    /** @brief Apply beta D beta-dagger, preserving pool collectives on empty ranks. */
    void apply(const UnitCell& cell, const pseudopot_cell_vnl& pp,
               const ModulePW::PW_Basis_K& basis, const int spin,
               const int bands, const int nbasis, const int npol, const int npw,
               const bool is_first_node, const T* vkb, const T* psi, T* hpsi);

    /** @brief Borrow the most recently computed projection coefficients. */
    T* get_becp() const { return becp_; }
};
} // namespace hamilt
#endif

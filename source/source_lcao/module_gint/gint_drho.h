#pragma once

#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

// Gint_drho integrates, on the real-space grid, the gradient density
//   [grad rho]^d(r) = sum_{Kk,Ll} D_{Kk,Ll} (grad^d phi_Kk)(r) phi_Ll(r),   d = x,y,z
// i.e. the derivative is taken on the FIRST (row) orbital of the density matrix.
// Pass a symmetrized matrix (D + D^T) to obtain
//   [grad rho^S]^d(r) = sum_{Kk,Ll} D_{Kk,Ll} [grad phi_Kk phi_Ll + phi_Kk grad phi_Ll].
//
// The grid loop / dm-to-gint preparation mirror Gint_rho; the orbital gradient
// (set_phi_dphi) is obtained exactly as in Gint_dvlocal. Everything is fp64 because
// set_phi_dphi only provides double-precision gradients.
//
// The three output buffers are ACCUMULATED into (phi_dot_phi uses +=), so the caller
// must zero-initialize drho_{x,y,z}[is] (each of length nrxx) before calling cal_gint().
class Gint_drho : public Gint
{
    public:
    Gint_drho(
        const std::vector<HContainer<double>*>& dm_vec,
        const int nspin,
        double** drho_x,
        double** drho_y,
        double** drho_z)
        : dm_vec_(dm_vec), nspin_(nspin),
          drho_x_(drho_x), drho_y_(drho_y), drho_z_(drho_z) {}

    void cal_gint();

    private:
    std::vector<HContainer<double>> init_dm_gint_() const;

    void cal_drho_(const std::vector<HContainer<double>>& dm_gint_vec) const;

    // input
    const std::vector<HContainer<double>*> dm_vec_;
    const int nspin_;

    // output: [grad rho]_{x,y,z}[is](ir), accumulated
    double** drho_x_ = nullptr;
    double** drho_y_ = nullptr;
    double** drho_z_ = nullptr;
};

}

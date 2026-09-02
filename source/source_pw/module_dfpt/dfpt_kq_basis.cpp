// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in design phase and has not been
// put into production yet.
// It may change in the future.
// Please use this code with caution.
// Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_kq_basis.h"

#include "source_base/global_function.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include <set>

namespace ModuleDFPT {

DFPT_KQ_Basis::DFPT_KQ_Basis() {}

DFPT_KQ_Basis::~DFPT_KQ_Basis() {}

void DFPT_KQ_Basis::init(const ModulePW::PW_Basis_K* pw_wfc,
                         const ModulePW::PW_Basis* pw_rho,
                         const ModuleBase::Vector3<double>& q_cart,
                         int ik)
{
    pw_wfc_ = pw_wfc;
    npwk_ = 0;
    ig_rho_.clear();
    gk2_.clear();
    gcar_.clear();

    if (pw_wfc_ == nullptr || pw_rho == nullptr)
    {
        return;
    }

    // DFPT couples k and k+q symmetrically; the perturbation wavevector is
    // generally incommensurate with the gamma-ladder, so the ground-state
    // basis must be a full complex basis (see class documentation).
    if (pw_wfc_->gamma_only)
    {
        ModuleBase::WARNING_QUIT("DFPT_KQ_Basis",
                                 "DFPT requires a complex (gamma_only=false) wavefunction basis for k+q. "
                                 "Please disable gamma_only for the wavefunction basis used by DFPT.");
    }

    // the two bases exchange G vectors through the shared FFT cell position
    if (pw_wfc_->nx != pw_rho->nx || pw_wfc_->ny != pw_rho->ny
        || pw_wfc_->nz != pw_rho->nz)
    {
        ModuleBase::WARNING_QUIT("DFPT_KQ_Basis",
                                 "DFPT requires the wavefunction and charge FFT grids to share "
                                 "their dimensions.");
    }

    const ModuleBase::Vector3<double> k_c = pw_wfc_->kvec_c[ik];
    kplusq_c_ = k_c + q_cart;

    // rho-grid reverse map of the shared FFT cell, used to attach the
    // charge-basis index of every enumerated k+q plane wave
    std::vector<int> ig_of_cell(pw_rho->nxyz, -1);
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        const int isz = pw_rho->ig2isz[ig];
        const int iz = isz % pw_rho->nz;
        const int is = isz / pw_rho->nz;
        const int ixy = pw_rho->is2fftixy[is];
        const int ix = ixy / pw_rho->fftny;
        const int iy = ixy % pw_rho->fftny;
        ig_of_cell[(ix * pw_rho->ny + iy) * pw_rho->nz + iz] = ig;
    }

    std::set<int> taken;
    auto try_push = [&](const int ix_in, const int iy_in, const int iz_in,
                        const ModuleBase::Matrix3& gbase,
                        const int ig_rho_hint) {
        int ix = ix_in;
        int iy = iy_in;
        int iz = iz_in;
        if (ix >= static_cast<int>(pw_wfc_->nx / 2) + 1)
        {
            ix -= pw_wfc_->nx;
        }
        if (iy >= static_cast<int>(pw_wfc_->ny / 2) + 1)
        {
            iy -= pw_wfc_->ny;
        }
        if (iz >= static_cast<int>(pw_wfc_->nz / 2) + 1)
        {
            iz -= pw_wfc_->nz;
        }
        const ModuleBase::Vector3<double> gcar
            = ModuleBase::Vector3<double>(ix, iy, iz) * gbase;
        const ModuleBase::Vector3<double> gpluskq = gcar + kplusq_c_;
        const double gk2 = gpluskq * gpluskq;
        if (gk2 > pw_wfc_->gk_ecut)
        {
            return;
        }
        const int cell = ((ix_in * pw_wfc_->ny) + iy_in) * pw_wfc_->nz + iz_in;
        if (!taken.insert(cell).second)
        {
            return;
        }
        gcar_.push_back(gcar);
        gk2_.push_back(gk2);
        ig_rho_.push_back(ig_rho_hint >= 0 ? ig_rho_hint : ig_of_cell[cell]);
    };

    // first reservoir: the wavefunction G grid (covers the ground-state
    // k-mesh balls)
    for (int ig = 0; ig < pw_wfc_->npw; ++ig)
    {
        const int isz = pw_wfc_->ig2isz[ig];
        const int iz = isz % pw_wfc_->nz;
        const int is = isz / pw_wfc_->nz;
        const int ixy = pw_wfc_->is2fftixy[is];
        const int ix = ixy / pw_wfc_->fftny;
        const int iy = ixy % pw_wfc_->fftny;
        try_push(ix, iy, iz, pw_wfc_->G, -1);
    }
    // completing reservoir: the charge G grid; its ball radius
    // (2*sqrt(ecutwfc) by construction) covers every q-shifted ball for q
    // inside the first Brillouin zone, which the wavefunction grid does not
    for (int ig = 0; ig < pw_rho->npw; ++ig)
    {
        const int isz = pw_rho->ig2isz[ig];
        const int iz = isz % pw_rho->nz;
        const int is = isz / pw_rho->nz;
        const int ixy = pw_rho->is2fftixy[is];
        const int ix = ixy / pw_rho->fftny;
        const int iy = ixy % pw_rho->fftny;
        try_push(ix, iy, iz, pw_rho->G, ig);
    }
    npwk_ = static_cast<int>(gcar_.size());
}

void DFPT_KQ_Basis::clear()
{
    pw_wfc_ = nullptr;
    kplusq_c_ = ModuleBase::Vector3<double>();
    npwk_ = 0;
    ig_rho_.clear();
    gk2_.clear();
    gcar_.clear();
}

} // namespace ModuleDFPT

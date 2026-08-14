// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_kq_basis.h"
#include "source_base/global_function.h"
#include "source_basis/module_pw/pw_basis_k.h"

namespace ModuleDFPT {

DFPT_KQ_Basis::DFPT_KQ_Basis() {}
DFPT_KQ_Basis::~DFPT_KQ_Basis() {}

void DFPT_KQ_Basis::init(const ModulePW::PW_Basis_K* pw_wfc,
                         const ModuleBase::Vector3<double>& q_cart,
                         int ik)
{
    pw_wfc_ = pw_wfc;
    npwk_ = 0;
    igl2ig_.clear();
    gk2_.clear();
    gcar_.clear();

    if (pw_wfc_ == nullptr)
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

    const ModuleBase::Vector3<double> k_c = pw_wfc_->kvec_c[ik];
    kplusq_c_ = k_c + q_cart;

    // Reuse the ground-state k-basis G grid (shared by all k points of the
    // pool): the k+q ball is a subset of it for every ik and q (see the
    // class documentation), so only the shifted-center selection is needed.
    const int npw = pw_wfc_->npw;
    for (int ig = 0; ig < npw; ++ig)
    {
        const int isz = pw_wfc_->ig2isz[ig];
        int iz = isz % pw_wfc_->nz;
        const int is = isz / pw_wfc_->nz;
        const int ixy = pw_wfc_->is2fftixy[is];
        int ix = ixy / pw_wfc_->fftny;
        int iy = ixy % pw_wfc_->fftny;
        if (ix >= int(pw_wfc_->nx / 2) + 1)
        {
            ix -= pw_wfc_->nx;
        }
        if (iy >= int(pw_wfc_->ny / 2) + 1)
        {
            iy -= pw_wfc_->ny;
        }
        if (iz >= int(pw_wfc_->nz / 2) + 1)
        {
            iz -= pw_wfc_->nz;
        }
        ModuleBase::Vector3<double> f(ix, iy, iz);
        const ModuleBase::Vector3<double> gcar = f * pw_wfc_->G;
        const ModuleBase::Vector3<double> gpluskq = gcar + kplusq_c_;
        const double gk2 = gpluskq * gpluskq;
        if (gk2 <= pw_wfc_->gk_ecut)
        {
            igl2ig_.push_back(ig);
            gcar_.push_back(gcar);
            gk2_.push_back(gk2);
        }
    }
    npwk_ = static_cast<int>(igl2ig_.size());
}

void DFPT_KQ_Basis::clear()
{
    pw_wfc_ = nullptr;
    kplusq_c_ = ModuleBase::Vector3<double>();
    npwk_ = 0;
    igl2ig_.clear();
    gk2_.clear();
    gcar_.clear();
}

int DFPT_KQ_Basis::get_ig2isz(int igl) const
{
    return pw_wfc_->ig2isz[igl2ig_[igl]];
}

} // namespace ModuleDFPT
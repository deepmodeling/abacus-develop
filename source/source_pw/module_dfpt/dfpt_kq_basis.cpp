#include "dfpt_kq_basis.h"

#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"

#include <string>

namespace ModuleDFPT
{

DFPT_KQ_Basis::DFPT_KQ_Basis()
{
}

DFPT_KQ_Basis::~DFPT_KQ_Basis()
{
}

void DFPT_KQ_Basis::init(const ModulePW::PW_Basis_K* pw_wfc,
                         const ModulePW::PW_Basis* pw_rho,
                         const ModuleBase::Vector3<double>& q_cart,
                         int ik)
{
    ModuleBase::TITLE("DFPT_KQ_Basis", "init");
    ModuleBase::timer::start("DFPT_KQ_Basis", "init");
    pw_wfc_ = pw_wfc;
    npwk_ = 0;
    ig_rho_.clear();
    gk2_.clear();
    gcar_.clear();

    if (pw_wfc_ == nullptr || pw_rho == nullptr)
    {
        ModuleBase::timer::end("DFPT_KQ_Basis", "init");
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
    if (pw_wfc_->nx != pw_rho->nx || pw_wfc_->ny != pw_rho->ny || pw_wfc_->nz != pw_rho->nz)
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

    // single-source enumeration over the full shared FFT cell: every integer
    // G with |G + (k+q)|^2 <= gk_ecut is a member of the k+q ball, and each
    // is carried by exactly one FFT cell (the half-space wrap makes the
    // cell-to-G map bijective), so scanning the whole cell space counts the
    // ball exactly -- no dependence on the ball layout of either
    // ground-state grid. A member whose cell carries no local charge-grid G
    // would break the veff convolution (kq2rho = -1); with ecutrho >= 4
    // ecutwfc and q inside the first Brillouin zone the charge ball has
    // twice the radius and carries every member, and the orphan counter
    // turns any violation into an explicit quit instead of a silent
    // coefficient zeroing.
    int orphan = 0;
    const double ball_cut2 = pw_wfc_->gk_ecut;
    const int hw_x = static_cast<int>(pw_wfc_->nx / 2) + 1;
    const int hw_y = static_cast<int>(pw_wfc_->ny / 2) + 1;
    const int hw_z = static_cast<int>(pw_wfc_->nz / 2) + 1;
    for (int ix = 0; ix < pw_wfc_->nx; ++ix)
    {
        const int wx = (ix >= hw_x) ? ix - pw_wfc_->nx : ix;
        for (int iy = 0; iy < pw_wfc_->ny; ++iy)
        {
            const int wy = (iy >= hw_y) ? iy - pw_wfc_->ny : iy;
            for (int iz = 0; iz < pw_wfc_->nz; ++iz)
            {
                const int wz = (iz >= hw_z) ? iz - pw_wfc_->nz : iz;
                const ModuleBase::Vector3<double> gcar = ModuleBase::Vector3<double>(wx, wy, wz) * pw_rho->G;
                const ModuleBase::Vector3<double> gpluskq = gcar + kplusq_c_;
                const double gk2 = gpluskq * gpluskq;
                if (gk2 > ball_cut2)
                {
                    continue;
                }
                const int ig = ig_of_cell[(ix * pw_wfc_->ny + iy) * pw_wfc_->nz + iz];
                if (ig < 0)
                {
                    ++orphan;
                    continue;
                }
                gcar_.push_back(gcar);
                gk2_.push_back(gk2);
                ig_rho_.push_back(ig);
            }
        }
    }
    if (orphan > 0)
    {
        ModuleBase::WARNING_QUIT("DFPT_KQ_Basis::init",
                                 "the k+q ball contains " + std::to_string(orphan)
                                     + " plane waves whose FFT cell carries no charge-grid G; "
                                       "increase ecutrho (>= 4 x ecutwfc covers any q inside the "
                                       "first Brillouin zone).");
    }
    npwk_ = static_cast<int>(gcar_.size());
    ModuleBase::timer::end("DFPT_KQ_Basis", "init");
}

void DFPT_KQ_Basis::clear()
{
    ModuleBase::TITLE("DFPT_KQ_Basis", "clear");
    ModuleBase::timer::start("DFPT_KQ_Basis", "clear");
    pw_wfc_ = nullptr;
    kplusq_c_ = ModuleBase::Vector3<double>();
    npwk_ = 0;
    ig_rho_.clear();
    gk2_.clear();
    gcar_.clear();
    ModuleBase::timer::end("DFPT_KQ_Basis", "clear");
}

} // namespace ModuleDFPT

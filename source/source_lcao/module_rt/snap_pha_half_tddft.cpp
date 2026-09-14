#include "snap_pha_half_tddft.h"

#include "source_base/vector3.h"
#include "source_basis/module_ao/orb_read.h"
#include "source_lcao/module_rt/snap_proj_half_tddft.h"

#include <complex>
#include <vector>

namespace module_rt
{

void snap_phialpha_half_tddft(const LCAO_Orbitals& orb,
                              std::vector<std::vector<std::complex<double>>>& nlm,
                              const ModuleBase::Vector3<double>& R1,
                              const int& T1,
                              const int& L1,
                              const int& m1,
                              const int& N1,
                              const ModuleBase::Vector3<double>& R0,
                              const ModuleBase::Vector3<double>& A,
                              const bool& calc_r)
{
    SnapIntegrationOptions options;
    snap_phialpha_half_tddft(orb, nlm, R1, T1, L1, m1, N1, R0, A, calc_r, options);
}

void snap_phialpha_half_tddft(const LCAO_Orbitals& orb,
                              std::vector<std::vector<std::complex<double>>>& nlm,
                              const ModuleBase::Vector3<double>& R1,
                              const int& T1,
                              const int& L1,
                              const int& m1,
                              const int& N1,
                              const ModuleBase::Vector3<double>& R0,
                              const ModuleBase::Vector3<double>& A,
                              const bool& calc_r,
                              const SnapIntegrationOptions& options)
{
    std::vector<ProjectorChannel> channels;
    const int lmax_alpha = orb.Alpha[0].getLmax();
    for (int L = 0; L <= lmax_alpha; ++L)
    {
        const int nchi_L = orb.Alpha[0].getNchi(L);
        for (int N = 0; N < nchi_L; ++N)
        {
            const auto& alpha_ln = orb.Alpha[0].PhiLN(L, N);
            ProjectorChannel channel;
            channel.l = L;
            channel.mesh = alpha_ln.getNr();
            channel.rcut = alpha_ln.getRcut();
            channel.radial_times_r = alpha_ln.getPsi_r();
            channel.radial_grid = alpha_ln.getRadial();
            channels.push_back(channel);
        }
    }

    snap_projector_half_tddft(orb, channels, nlm, R1, T1, L1, m1, N1, R0, A, calc_r, options, "snap_phialpha_half_tddft");
}

} // namespace module_rt

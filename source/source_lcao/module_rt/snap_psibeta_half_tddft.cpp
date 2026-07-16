#include "snap_psibeta_half_tddft.h"

namespace module_rt
{

void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0,
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r)
{
    SnapIntegrationOptions options;
    snap_psibeta_half_tddft(orb, infoNL_, nlm, R1, T1, L1, m1, N1, R0, T0, A, calc_r, options);
}

void snap_psibeta_half_tddft(const LCAO_Orbitals& orb,
                             const InfoNonlocal& infoNL_,
                             std::vector<std::vector<std::complex<double>>>& nlm,
                             const ModuleBase::Vector3<double>& R1,
                             const int& T1,
                             const int& L1,
                             const int& m1,
                             const int& N1,
                             const ModuleBase::Vector3<double>& R0,
                             const int& T0,
                             const ModuleBase::Vector3<double>& A,
                             const bool& calc_r,
                             const SnapIntegrationOptions& options)
{
    std::vector<ProjectorChannel> channels;
    channels.reserve(infoNL_.nproj[T0]);

    // UPF nonlocal beta projectors already follow the r * beta_l(r) convention.
    for (int ip = 0; ip < infoNL_.nproj[T0]; ++ip)
    {
        const auto& proj = infoNL_.Beta[T0].Proj[ip];
        ProjectorChannel channel;
        channel.l = proj.getL();
        channel.mesh = proj.getNr();
        channel.dk = proj.getDk();
        channel.rcut = proj.getRcut();
        channel.radial_times_r = proj.getBeta_r();
        channel.radial_grid = proj.getRadial();
        channels.push_back(channel);
    }

    snap_projector_half_tddft(orb, channels, nlm, R1, T1, L1, m1, N1, R0, A, calc_r, options, "snap_psibeta_half_tddft");
}

} // namespace module_rt

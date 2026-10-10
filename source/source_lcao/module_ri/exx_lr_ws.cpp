#include "exx_lr_ws.h"
#include "exx_lri.h"

namespace
{
template<typename Tdata>
void copy_pack(const RI::Exx<int, int, 3, Tdata>& source,
               RI::Exx<int, int, 3, Tdata>& destination,
               const std::string& name)
{
    destination.lri.data_pool.emplace(name, source.lri.data_pool.at(name));
}
}

template<typename Tdata>
void share_exx_geometry(const RI::Exx<int, int, 3, Tdata>& source,
                       RI::Exx<int, int, 3, Tdata>& destination)
{
    // Only geometry packs may cross the KS/LR boundary. In particular, do not copy
    // Ds, Hs, cvc, label bindings, symmetry filters or parallel/post-processing objects.
    if (source.flag_finish.Cs) { copy_pack(source, destination, "Cs_"); }
    if (source.flag_finish.Vs) { copy_pack(source, destination, "Vs_"); }
    for (int axis = 0; axis < 3; ++axis)
    {
        const std::string suffix = std::to_string(axis) + "_";
        const std::string dc_name = "dCs_" + suffix;
        const std::string dv_name = "dVs_" + suffix;
        if (source.flag_finish.dCs) { copy_pack(source, destination, dc_name); }
        if (source.flag_finish.dVs) { copy_pack(source, destination, dv_name); }
        for (int second = 0; second < 3; ++second)
        {
            const std::string stress_suffix = suffix + std::to_string(second) + "_";
            const std::string dcr_name = "dCRs_" + stress_suffix;
            const std::string dvr_name = "dVRs_" + stress_suffix;
            if (source.flag_finish.dCRs) { copy_pack(source, destination, dcr_name); }
            if (source.flag_finish.dVRs) { copy_pack(source, destination, dvr_name); }
        }
    }
    destination.flag_finish.Cs = source.flag_finish.Cs;
    destination.flag_finish.Vs = source.flag_finish.Vs;
    destination.flag_finish.dCs = source.flag_finish.dCs;
    destination.flag_finish.dVs = source.flag_finish.dVs;
    destination.flag_finish.dCRs = source.flag_finish.dCRs;
    destination.flag_finish.dVRs = source.flag_finish.dVRs;
}

template<typename Tdata>
std::shared_ptr<Exx_LRI<Tdata>> Exx_LRI<Tdata>::make_lr_workspace(const UnitCell& ucell,
                                                              const K_Vectors& kv) const
{
    auto workspace = std::make_shared<Exx_LRI<Tdata>>(this->info);
    workspace->mpi_comm = this->mpi_comm;
    workspace->p_kv = &kv;
    workspace->abfs_Lmax_ = this->abfs_Lmax_;
    std::map<TA, TatomR> atoms_pos;
    for (int atom = 0; atom < ucell.nat; ++atom)
    {
        const int type = ucell.iat2it[atom];
        const int index = ucell.iat2ia[atom];
        atoms_pos[atom] = RI_Util::Vector3_to_array3(ucell.atoms[type].tau[index]);
    }
    const std::array<TatomR, Ndim> lattice = {RI_Util::Vector3_to_array3(ucell.a1),
                                            RI_Util::Vector3_to_array3(ucell.a2),
                                            RI_Util::Vector3_to_array3(ucell.a3)};
    const std::array<Tcell, Ndim> period = {kv.nmp[0], kv.nmp[1], kv.nmp[2]};
    workspace->exx_lri.set_parallel(this->mpi_comm, atoms_pos, lattice, period);
    share_exx_geometry(this->exx_lri, workspace->exx_lri);
    return workspace;
}

template void share_exx_geometry(const RI::Exx<int, int, 3, double>&,
                                RI::Exx<int, int, 3, double>&);
template void share_exx_geometry(const RI::Exx<int, int, 3, std::complex<double>>&,
                                RI::Exx<int, int, 3, std::complex<double>>&);
template std::shared_ptr<Exx_LRI<double>> Exx_LRI<double>::make_lr_workspace(const UnitCell&,
                                                                         const K_Vectors&) const;
template std::shared_ptr<Exx_LRI<std::complex<double>>>
Exx_LRI<std::complex<double>>::make_lr_workspace(const UnitCell&, const K_Vectors&) const;

#include "vxc_op_r.h"

#include "source_io/module_hs/hs_sparse_io.h"
#include "source_base/module_out/filename.h"

namespace ModuleIO
{

/// @brief Helper to calculate sparse HR representation
template <typename T>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, T>>> cal_HR_sparse(
    const hamilt::HContainer<T>& hR,
    const double sparse_thr)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, T>>> target;
    sparse_format::cal_HContainer<T>(*hR.get_paraV(), sparse_thr, hR, target);
    return target;
}

template <typename TK, typename TR>
void write_Vxc_R(const int nspin,
                 const Parallel_Orbitals* pv,
                 const UnitCell& ucell,
                 Structure_Factor& sf,
                 surchem& solvent,
                 const ModulePW::PW_Basis& rho_basis,
                 const ModulePW::PW_Basis& rhod_basis,
                 const ModuleBase::matrix& vloc,
                 const Charge& chg,
                 const K_Vectors& kv,
                 const std::vector<double>& orb_cutoff,
                 Grid_Driver& gd,
                 const std::string& global_out_dir,
                 bool cal_exx,
                 double hybrid_alpha,
                 bool real_number
#ifdef __EXX
                 ,
                 const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
                 const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
                 ,
                 const double sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "write_Vxc_R");

    // 1. real-space xc potential
    double etxc = 0.0;
    double vtxc = 0.0;
    elecstate::Potential potxc(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc);
    std::vector<std::string> compnents_list = {"xc"};
    potxc.pot_register(compnents_list);
    potxc.update_from_charge(&chg, &ucell);

    // 2. allocate H(R)
    // (the number of hR: 1 for nspin=1, 4; 2 for nspin=2)
    int nspin0 = (nspin == 2) ? 2 : 1;
    std::vector<hamilt::HContainer<TR>> vxcs_R_ao(nspin0, hamilt::HContainer<TR>(ucell, pv));
#ifdef __EXX
    std::array<int, 3> Rs_period = {kv.nmp[0], kv.nmp[1], kv.nmp[2]};
    const auto cell_nearest = hamilt::init_cell_nearest(ucell, Rs_period);
#endif
    for (int is = 0; is < nspin0; ++is)
    {
        if (std::is_same<TK, double>::value)
        {
            vxcs_R_ao[is].fix_gamma();
        }
#ifdef __EXX
        if (cal_exx)
        {
            real_number
                ? hamilt::reallocate_hcontainer(*Hexxd, &vxcs_R_ao[is], &cell_nearest)
                : hamilt::reallocate_hcontainer(*Hexxc, &vxcs_R_ao[is], &cell_nearest);
        }
#endif
    }

    // 3. calculate the Vxc(R)
    hamilt::HS_Matrix_K<TK> vxc_k_ao(pv, 1); // only hk is needed, sk is skipped
    for (int is = 0; is < nspin0; ++is)
    {
        hamilt::Veff<hamilt::OperatorLCAO<TK, TR>> vxcs_op_ao(&vxc_k_ao,
                                                              kv.kvec_d,
                                                              &potxc,
                                                              &vxcs_R_ao[is],
                                                              &ucell,
                                                              orb_cutoff,
                                                              &gd,
                                                              nspin);
        vxcs_op_ao.set_current_spin(is);
        vxcs_op_ao.contributeHR();
#ifdef __EXX
        if (cal_exx)
        {
            real_number ? RI_2D_Comm::add_HexxR(is,
                                                 hybrid_alpha,
                                                 *Hexxd,
                                                 *pv,
                                                 ucell.get_npol(),
                                                 vxcs_R_ao[is],
                                                 &cell_nearest)
                        : RI_2D_Comm::add_HexxR(is,
                                                 hybrid_alpha,
                                                 *Hexxc,
                                                 *pv,
                                                 ucell.get_npol(),
                                                 vxcs_R_ao[is],
                                                 &cell_nearest);
        }
#endif
    }

    // 4. write Vxc(R) in csr format
    for (int is = 0; is < nspin0; ++is)
    {
        std::set<Abfs::Vector3_Order<int>> all_R_coor = sparse_format::get_R_range(vxcs_R_ao[is]);
        const std::string filename = "Vxc_R_spin" + std::to_string(is);
        ModuleIO::SparseWriteOptions options;
        options.filename = global_out_dir + filename + ".csr";
        options.label = filename;
        options.threshold = sparse_thr;
        options.binary = false;
        options.istep = -1;
        options.reduce = true;
        ModuleIO::save_sparse(cal_HR_sparse(vxcs_R_ao[is], sparse_thr),
                              all_R_coor,
                              *pv,
                              options);
    }
}

// Explicit template instantiations
template void write_Vxc_R<double, double>(
    const int, const Parallel_Orbitals*, const UnitCell&, Structure_Factor&, surchem&,
    const ModulePW::PW_Basis&, const ModulePW::PW_Basis&, const ModuleBase::matrix&,
    const Charge&, const K_Vectors&, const std::vector<double>&, Grid_Driver&,
    const std::string&, bool, double, bool
#ifdef __EXX
    , const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
    , const double);

template void write_Vxc_R<std::complex<double>, double>(
    const int, const Parallel_Orbitals*, const UnitCell&, Structure_Factor&, surchem&,
    const ModulePW::PW_Basis&, const ModulePW::PW_Basis&, const ModuleBase::matrix&,
    const Charge&, const K_Vectors&, const std::vector<double>&, Grid_Driver&,
    const std::string&, bool, double, bool
#ifdef __EXX
    , const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
    , const double);

template void write_Vxc_R<std::complex<double>, std::complex<double>>(
    const int, const Parallel_Orbitals*, const UnitCell&, Structure_Factor&, surchem&,
    const ModulePW::PW_Basis&, const ModulePW::PW_Basis&, const ModuleBase::matrix&,
    const Charge&, const K_Vectors&, const std::vector<double>&, Grid_Driver&,
    const std::string&, bool, double, bool
#ifdef __EXX
    , const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>*,
    const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>*
#endif
    , const double);

} // namespace ModuleIO

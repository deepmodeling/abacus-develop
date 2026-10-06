#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa)
{
    ModuleBase::TITLE("RPA_LRI", "output_cut_coulomb_cs");
    ModuleBase::timer::start("RPA_LRI", "output_cut_coulomb_cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    std::cout << "Use rpa_ccp_rmesh_times=" << this->ccp_rmesh_times_cut << " to calculate cut Coulomb" << std::endl;
    // Shrink_ABFS_ORBITAL cannot exceed this angular momentum of MGT
    exx_lri_rpa->cal_cut_coulomb_cs(Vs_cut_IJR, Cs, ucell, this->runtime.input.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_cut, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJ
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_cut_IJR, atoms00, atoms01);
    Vs_cut_IJR.clear();
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJ, period);
    if (this->runtime.input.out_librpa_ver == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_lri_rpa,
                                  use_shrink ? "aux_basis_s.txt" : "aux_basis.txt",
                                  use_shrink ? "wfc_basis_s.txt" : "basis_map.txt");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "V_cut_", exx_lri_rpa);
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_cut_", exx_lri_rpa);
    }
    Vs_period.clear();
    Vs_period.swap(tmp);

    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_lri_rpa->exx_lri.post_2D.set_tensors_map2(this->Cs_period);

    if (this->runtime.input.out_librpa_ver == 1)
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
        {
            this->out_Cs_v1(ucell, this->Cs_period, "Cs_shrink_");
        }
        else
        {
            this->out_Cs_v1(ucell, this->Cs_period, "Cs_");
        }
    }
    else
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
            this->out_Cs(ucell, this->Cs_period, "Cs_shrinked_data_");
        else
            this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    Cs_period.clear();
    Cs_period.swap(tmp);

    ModuleBase::timer::end("RPA_LRI", "output_cut_coulomb_cs");
}

template <typename T, typename Tdata>
Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<T, Tdata>::select_coulomb_basis_method_(Exx_LRI<double>* exx_lri) const
{
    if (exx_lri == nullptr || exx_lri->exx_objs.empty())
    {
        throw std::invalid_argument("Cannot select Coulomb basis method from an empty Exx_LRI object.");
    }
    return exx_lri->exx_objs.count(Conv_Coulomb_Pot_K::Coulomb_Method::Center2)
               ? Conv_Coulomb_Pot_K::Coulomb_Method::Center2
               : exx_lri->exx_objs.begin()->first;
}

template <typename T, typename Tdata>
std::vector<int> RPA_LRI<T, Tdata>::collect_atom_naux_(const UnitCell& ucell, Exx_LRI<double>* exx_lri) const
{
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> atom_naux(static_cast<std::size_t>(ucell.nat), 0);
    for (int I = 0; I != ucell.nat; ++I)
    {
        atom_naux[static_cast<std::size_t>(I)]
            = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
        if (atom_naux[static_cast<std::size_t>(I)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 output found a non-positive per-atom auxiliary size.");
        }
    }
    return atom_naux;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "output_ewald_coulomb");
    ModuleBase::timer::start("RPA_LRI", "output_ewald_coulomb");

    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = this->ccp_rmesh_times_ewald;
    if (!exx_full_coulomb)
    {
        Exx_LRI<double>* new_exx = new Exx_LRI<double>(local_info);
        exx_full_coulomb.reset(new_exx);
    }

    if (this->info.shrink_abfs_pca_thr >= 0.0)
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs_shrink);
    else
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    exx_full_coulomb->cal_ewald_coulomb(Vs_full_IJR, Cs, ucell, this->runtime.input.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_ewald, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, Vs_full_IJR, atoms00, atoms01);
    Vs_full_IJR.clear();

    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_full_IJ, period);
    if (this->runtime.input.out_librpa_ver == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_full_coulomb.get(),
                                  use_shrink ? "aux_basis_s.txt" : "aux_basis.txt",
                                  use_shrink ? "wfc_basis_s.txt" : "basis_map.txt");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "V_full_", exx_full_coulomb.get());
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_mat_", exx_full_coulomb.get());
    }
    Vs_period.clear();
    Vs_period.swap(tmp);
    Cs.clear();
    Cs.swap(tmp);

    exx_full_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "output_ewald_coulomb");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k(const UnitCell& ucell,
                                      std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                      std::string filename,
                                      Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k");
    ModuleBase::timer::start("RPA_LRI", "out_coulomb_k");

    int all_mu = 0;
    std::vector<int> mu_shift(ucell.nat);
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (this->runtime.rank + 1) << ".txt";

    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Vs)
    {
        auto I = Ip.first;
        size_t mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;

                auto R = JPp.first.second;
                if (J < I)
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                auto iJ = vq_Jp.first;
                auto& vq_J = vq_Jp.second;
                size_t nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[iJ]);
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[iJ] + 1
                    << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * this->runtime.input.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    ofs << std::setw(25) << (*vq_J.data)[i].real() << std::setw(25) << (*vq_J.data)[i].imag()
                        << std::endl;
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_coulomb_k");
}

template void RPA_LRI<double, double>::output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa);
template void RPA_LRI<std::complex<double>, double>::output_cut_coulomb_cs(const UnitCell& ucell,
                                                                           Exx_LRI<double>* exx_lri_rpa);
template Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<double, double>::select_coulomb_basis_method_(
    Exx_LRI<double>* exx_lri) const;
template Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<std::complex<double>, double>::select_coulomb_basis_method_(
    Exx_LRI<double>* exx_lri) const;
template std::vector<int> RPA_LRI<double, double>::collect_atom_naux_(const UnitCell& ucell,
                                                                      Exx_LRI<double>* exx_lri) const;
template std::vector<int> RPA_LRI<std::complex<double>, double>::collect_atom_naux_(const UnitCell& ucell,
                                                                                    Exx_LRI<double>* exx_lri) const;
template void RPA_LRI<double, double>::output_ewald_coulomb(const UnitCell& ucell,
                                                            const K_Vectors& kv,
                                                            const LCAO_Orbitals& orb);
template void RPA_LRI<std::complex<double>, double>::output_ewald_coulomb(const UnitCell& ucell,
                                                                          const K_Vectors& kv,
                                                                          const LCAO_Orbitals& orb);
template void RPA_LRI<double, double>::out_coulomb_k(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Vs,
    std::string filename,
    Exx_LRI<double>* exx_lri);
template void RPA_LRI<std::complex<double>, double>::out_coulomb_k(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Vs,
    std::string filename,
    Exx_LRI<double>* exx_lri);

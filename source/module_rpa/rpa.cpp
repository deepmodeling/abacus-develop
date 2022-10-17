#include "rpa.h"

#include "../src_ri/conv_coulomb_pot-inl.h"
#include "../src_ri/conv_coulomb_pot_k-template.h"
#include "../src_ri/exx_abfs-abfs_index.h"
#include "../src_ri/exx_abfs-construct_orbs.h"
#include "../src_ri/exx_abfs-io.h"
#include "../src_ri/exx_abfs-util.h"
#include "src_lcao/global_fp.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <mpi.h>

namespace ModuleRPA
{
void DFT_RPA_interface::out_for_RPA(const Parallel_Orbitals &parav,
                                    const psi::Psi<std::complex<double>> &psi,
                                    Local_Orbital_Charge &loc)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_for_RPA");
    this->out_bands();
    this->out_eigen_vector(parav, psi);
    this->out_struc();

    rpa_exx_lcao_.exx_init();
    std::cout << "rpa_pca_threshold: " << rpa_exx_lcao_.info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times: " << rpa_exx_lcao_.info.ccp_rmesh_times << std::endl;
    rpa_exx_lcao_.cal_exx_ions(parav);
    rpa_exx_lcao_.cal_exx_elec(loc, nullptr);
    
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << rpa_exx_lcao_.get_energy() / 2.0
              << std::endl;
    this->out_Cs();
    this->out_coulomb_k();

    std::cout << "etxc(Ha):" << std::fixed << std::setprecision(15) << GlobalC::en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha):" << std::fixed << std::setprecision(15) << GlobalC::en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha):" << std::fixed << std::setprecision(15)
              << (GlobalC::en.etot - GlobalC::en.etxc + rpa_exx_lcao_.get_energy()) / 2.0 << std::endl;

    return;
}

void DFT_RPA_interface::out_eigen_vector(const Parallel_Orbitals &parav, const psi::Psi<std::complex<double>> &psi)
{

    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = GlobalV::NSPIN == 2 ? GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    const int npsin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    const std::complex<double> zero(0.0, 0.0);

    for (int ik = 0; ik < nks_tot; ik++)
    {
        std::stringstream ss;
        ss << "KS_eigenvector_" << ik << ".dat";

        std::ofstream ofs;
        if (GlobalV::MY_RANK == 0)
            ofs.open(ss.str().c_str(), std::ios::out);
        std::vector<ModuleBase::ComplexMatrix> is_wfc_ib_iw(npsin_tmp);
        for (int is = 0; is < npsin_tmp; is++)
        {
            ofs << ik + 1 << std::endl;
            
            is_wfc_ib_iw[is].create(GlobalV::NBANDS,GlobalV::NLOCAL);
            for (int ib_global = 0; ib_global < GlobalV::NBANDS; ++ib_global)
            {
                std::vector<std::complex<double>> wfc_iks(GlobalV::NLOCAL, zero);

                const int ib_local = parav.trace_loc_col[ib_global];

                if (ib_local >= 0)
                    for (int ir = 0; ir < psi.get_nbasis(); ir++)
                        wfc_iks[parav.MatrixInfo.row_set[ir]] = psi(ik + nks_tot * is, ib_local, ir);

                std::vector<std::complex<double>> tmp = wfc_iks;
#ifdef __MPI
                MPI_Allreduce(&tmp[0], &wfc_iks[0], GlobalV::NLOCAL, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
                for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                    is_wfc_ib_iw[is](ib_global,iw)=wfc_iks[iw];
            } // ib
        } // is
        for(int iw=0; iw<GlobalV::NLOCAL; iw++)
        {
            for(int ib=0; ib<GlobalV::NBANDS; ib++)
            {
                for(int is=0; is < npsin_tmp; is++)
                    ofs << std::setw(21) << std::fixed << std::setprecision(15) << is_wfc_ib_iw[is](ib,iw).real() << std::setw(21)
                        << std::fixed << std::setprecision(15) << is_wfc_ib_iw[is](ib,iw).imag() << std::endl;
            }
        }
        ofs.close();
    } // ik
    return;
}

void DFT_RPA_interface::out_struc()
{
    if (GlobalV::MY_RANK != 0)
        return;
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    double TWOPI_Bohr2A = ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    ModuleBase::Matrix3 lat = GlobalC::ucell.latvec / ModuleBase::BOHR_TO_A;
    ModuleBase::Matrix3 G = GlobalC::ucell.G * TWOPI_Bohr2A;
    std::stringstream ss;
    ss << "stru_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G.e11 << std::setw(15) << G.e12 << std::setw(15) << G.e13 << std::endl;
    ofs << G.e21 << std::setw(15) << G.e22 << std::setw(15) << G.e23 << std::endl;
    ofs << G.e31 << std::setw(15) << G.e32 << std::setw(15) << G.e33 << std::endl;

    ofs << GlobalC::kv.nmp[0] << std::setw(6) << GlobalC::kv.nmp[1] << std::setw(6) << GlobalC::kv.nmp[2]
        << std::setw(6) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
        ofs << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].x * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].y * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].z * TWOPI_Bohr2A
            << std::endl;
    ofs.close();
    return;
}

void DFT_RPA_interface::out_bands()
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (GlobalV::MY_RANK != 0)
        return;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    const int nspin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    std::stringstream ss;
    ss << "band_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << GlobalC::kv.nks << std::endl;
    ofs << GlobalV::NSPIN << std::endl;
    ofs << GlobalV::NBANDS << std::endl;
    ofs << GlobalV::NLOCAL << std::endl;
    ofs << (GlobalC::en.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != GlobalV::NBANDS; ib++)
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << GlobalC::wf.wg(ik, ib) * nks_tot
                    << std::setw(18) << std::fixed << std::setprecision(8)
                    << GlobalC::wf.ekb[ik + is * nks_tot][ib] / 2.0 << std::setw(18) << std::fixed
                    << std::setprecision(8) << GlobalC::wf.ekb[ik + is * nks_tot][ib] * ModuleBase::Ry_to_eV
                    << std::endl;
        }
    }
    ofs.close();
    return;
}

void DFT_RPA_interface::out_Cs()
{
    auto cal_atom_centres_core = [](const std::vector<std::pair<size_t, size_t>> &atom_pairs_core) -> std::set<size_t> {
        std::set<size_t> atom_centres_core;
        for (const std::pair<size_t, size_t> &atom_pair: atom_pairs_core)
        {
            atom_centres_core.insert(atom_pair.first);
            atom_centres_core.insert(atom_pair.second);
        }
        return atom_centres_core;
    };

    const std::set<size_t> atom_centres_core = cal_atom_centres_core(rpa_exx_lcao_.get_atom_pairs_core());

    std::stringstream ss;
    ss << "Cs_data.txt";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << atom_centres_core.size() << "    " << 0 << std::endl;
    for (auto &Ip: rpa_exx_lcao_.get_Cps())
    {
        size_t I = Ip.first;
        size_t i_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[I]].nw;
        size_t mu_num = rpa_exx_lcao_.get_index_abfs()[GlobalC::ucell.iat2it[I]].count_size;
        for (auto &Jp: Ip.second)
        {
            size_t J = Jp.first;
            size_t j_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[J]].nw;
            for (auto &Rp: Jp.second)
            {
                auto R = Rp.first;
                ofs << I + 1 << "   " << J + 1 << "   " << R.x << "   " << R.y << "   " << R.z << "   " << i_num
                    << std::endl;
                ofs << j_num << "   " << mu_num << std::endl;
                auto &tmp_Cs = *Rp.second;
                size_t ele_num = tmp_Cs.nr * tmp_Cs.nc;
                for (int it = 0; it != ele_num; it++)
                    ofs << std::setw(15) << std::fixed << std::setprecision(9) << tmp_Cs.c[it] << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

void DFT_RPA_interface::out_coulomb_k()
{
    int all_mu = 0;
    vector<int> mu_shift(GlobalC::ucell.nat);
    for (int I=0;I!=GlobalC::ucell.nat;I++)
    {
        mu_shift[I] = all_mu;
        all_mu += rpa_exx_lcao_.get_index_abfs()[GlobalC::ucell.iat2it[I]].count_size;
    }

    std::stringstream ss;
    ss << "coulomb_mat_"<<GlobalV::MY_RANK<<".txt";

    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);

    ofs << GlobalC::kv.nks << std::endl;
    for (auto &Ip: rpa_exx_lcao_.get_Vps())
    {
        auto I = Ip.first;
        size_t mu_num = rpa_exx_lcao_.get_index_abfs()[GlobalC::ucell.iat2it[I]].count_size;
        for (auto &Jp: Ip.second)
        {
            auto J = Jp.first;
            size_t nu_num = rpa_exx_lcao_.get_index_abfs()[GlobalC::ucell.iat2it[J]].count_size;

            for (int ik = 0; ik != GlobalC::kv.nks; ik++)
            {
                ModuleBase::ComplexMatrix tmp_Vk(mu_num, nu_num);
                tmp_Vk.zero_out();
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[J] + 1
                    << "   " << mu_shift[J] + nu_num << std::endl;
                ofs << ik + 1 << "  " << GlobalC::kv.wk[ik] / 2.0 << std::endl;
                for (auto &Rp: Jp.second)
                {
                    auto R = Rp.first;
                    ModuleBase::ComplexMatrix tmp_VR(*Rp.second);
                    const double arg
                        = 1 * (GlobalC::kv.kvec_c[ik] * (R * GlobalC::ucell.latvec)) * ModuleBase::TWO_PI; // latvec
                    const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                    tmp_Vk += tmp_VR * kphase;
                }
                for (int i = 0; i != tmp_Vk.size; i++)
                {
                    ofs << std::setw(21) << std::fixed << std::setprecision(12) << tmp_Vk.c[i].real() << std::setw(21)
                        << std::fixed << std::setprecision(12) << tmp_Vk.c[i].imag() << std::endl;
                }
            }
        }
    }
    ofs.close();
}

void RPAExxLcao::exx_init()
{
    std::cout << "rpa_exx_init!!!" << std::endl;
#ifdef __MPI
    if (GlobalC::exx_global.info.separate_loop)
    {
        Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
        Hexx_para.mixing_beta = 0;
    }
    else
    {
        if ("plain" == GlobalC::CHR.mixing_mode)
        {
            Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
        }
        else if ("pulay" == GlobalC::CHR.mixing_mode)
        {
            Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
        }
        else
        {
            throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
        }
        Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
    }
#endif

    lcaos = Exx_Abfs::Construct_Orbs::change_orbs(GlobalC::ORB, kmesh_times);
#ifdef __MPI
    Exx_Abfs::Util::bcast(info.files_abfs, 0, MPI_COMM_WORLD);
#endif

    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_same_atom
        = Exx_Abfs::Construct_Orbs::abfs_same_atom(lcaos, kmesh_times, info.pca_threshold);
    if (info.files_abfs.empty())
    {
        abfs = abfs_same_atom;
    }
    else
    {
        abfs = Exx_Abfs::IO::construct_abfs(abfs_same_atom, GlobalC::ORB, info.files_abfs, kmesh_times);
    }

    switch (info.hybrid_type)
    {
    case Exx_Global::Hybrid_Type::HF:
        abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hf, {}, info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::No:
        abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hf, {}, info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::PBE0:
        abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hf, {}, info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::HSE:
        abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(abfs,
                                                    Conv_Coulomb_Pot_K::Ccp_Type::Hse,
                                                    {{"hse_omega", info.hse_omega}},
                                                    info.ccp_rmesh_times);
        break;
    default:
        throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
    }

    for (size_t T = 0; T != abfs.size(); ++T)
    {
        Exx_Abfs::Lmax = std::max(Exx_Abfs::Lmax, static_cast<int>(abfs[T].size()) - 1);
    }

    const ModuleBase::Element_Basis_Index::Range &&range_lcaos = Exx_Abfs::Abfs_Index::construct_range(lcaos);
    index_lcaos = ModuleBase::Element_Basis_Index::construct_index(range_lcaos);

    const ModuleBase::Element_Basis_Index::Range &&range_abfs = Exx_Abfs::Abfs_Index::construct_range(abfs);
    index_abfs = ModuleBase::Element_Basis_Index::construct_index(range_abfs);

    m_abfs_abfs.init(2, kmesh_times, (1 + info.ccp_rmesh_times) / 2.0);

    m_abfs_abfs.init_radial(abfs_ccp, abfs);

    m_abfslcaos_lcaos.init(1, kmesh_times, 1);

    m_abfslcaos_lcaos.init_radial(abfs_ccp, lcaos, lcaos);

    Born_von_Karman_period = ModuleBase::Vector3<int>{GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]};
}

/*
void RPAExxLcao::exx_cal_ions()
{
    std::cout<<"rpa_exx_cal_ions!!!"<<std::endl;
    auto cal_atom_centres_core = [](const std::vector<std::pair<size_t, size_t>> &atom_pairs_core) -> std::set<size_t> {
        std::set<size_t> atom_centres_core;
        for (const std::pair<size_t, size_t> &atom_pair: atom_pairs_core)
        {
            atom_centres_core.insert(atom_pair.first);
            atom_centres_core.insert(atom_pair.second);
        }
        return atom_centres_core;
    };

#ifdef __MPI
    if (atom_pairs_core_origin.empty()) {
        switch (info.distribute_type)
        {
        case Exx_Lcao::Distribute_Type::Htime:
            atom_pairs_core_origin
                = Exx_Abfs::Parallel::Distribute::Htime::distribute(Born_von_Karman_period, info.ccp_rmesh_times);
            break;
        case Exx_Lcao::Distribute_Type::Kmeans2:
            atom_pairs_core_origin
                = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans2(MPI_COMM_WORLD);
            break;
        case Exx_Lcao::Distribute_Type::Kmeans1:
            atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans1(
                MPI_COMM_WORLD, info.ccp_rmesh_times);
            break;
        case Exx_Lcao::Distribute_Type::Order:
            atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Order::distribute(info.ccp_rmesh_times);
            break;
        default:
            throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                    + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        }
    }
#endif

    init_radial_table_ions(cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin);

    Vs = Abfs::cal_Vs(atom_pairs_core_origin, m_abfs_abfs, index_abfs, info.ccp_rmesh_times, info.v_threshold, Vws);

    Abfs::delete_empty_ptrs(Vws);

    Vps = Abfs::cal_mps(Born_von_Karman_period, Vs);

    atom_pairs_core = Abfs::get_atom_pair(Vps);

    const std::set<size_t> atom_centres_core = cal_atom_centres_core(atom_pairs_core);

    H_atom_pairs_core = Abfs::get_H_pairs_core(atom_pairs_core);

    Abfs::delete_empty_ptrs(Cws);
}
*/

void DFT_RPA_interface::print_matrix(char *desc, const ModuleBase::matrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf(" %9.5f", mat.c[i * nc + j]);
        printf("\n");
    }
}

void DFT_RPA_interface::print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf("%10.6f,%8.6f ", mat.c[i * nc + j].real(), mat.c[i * nc + j].imag());
        printf("\n");
    }
}

} // namespace ModuleRPA

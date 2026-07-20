//=======================
// AUTHOR : jiyy
// DATE :   2024-03-08
//=======================

#ifndef EWALD_VQ_HPP
#define EWALD_VQ_HPP

#include <RI/comm/mix/Communicate_Tensors_Map_Judge.h>
#include <RI/distribute/Divide_Atoms.h>
#include <RI/distribute/Distribute_Equally.h>
#include <RI/global/Global_Func-1.h>

// #include <chrono>
#include "ri_2d_comm.h"
#include "ri_util.h"
#include "conv_coulomb_pot_k.h"
#include "ewald_mpi_utils.h"
#include "exx_abfs-construct_orbs.h"
#include "exx_rotate_abfs.h"
#include "gaussian_abfs.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_base/element_basis_index.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "singular_value.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace EwaldVqDetail
{
template<typename T>
inline RI::Tensor<T> tensor_adjoint(const RI::Tensor<T>& tensor)
{
    return tensor.transpose();
}

inline RI::Tensor<std::complex<double>> tensor_adjoint(const RI::Tensor<std::complex<double>>& tensor)
{
    return tensor.dagger();
}
} // namespace EwaldVqDetail

template<typename Tdata>
Ewald_Vq<Tdata>::Ewald_Vq(){}

template<typename Tdata>
Ewald_Vq<Tdata>::~Ewald_Vq(){}

template <typename Tdata>
void Ewald_Vq<Tdata>::init(const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const MPI_Comm& mpi_comm_in,
                           const K_Vectors* kv_in,
                           std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& lcaos_in,
                           std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in,
                           const std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string,std::string>>> &coulomb_param_in,
                           std::shared_ptr<ORB_gaunt_table> MGT_in,
                           const double &ccp_rmesh_times_in,
                           const double &exact_ccp_rmesh_times_in,
                           const double &ewald_lambda_in,
                           const double &kmesh_times_in,
                           const int& ewald_dimension_in,
                           const ModuleBase::Element_Basis_Index::IndexPermutation &abfs_old_to_new)
{
    ModuleBase::TITLE("Ewald_Vq", "init");
    ModuleBase::timer::start("Ewald_Vq", "init");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = kv_in;
    this->nks0 = this->p_kv->get_nkstot_nospin();
    this->kvec_c.resize(this->nks0);
    this->ccp_rmesh_times = ccp_rmesh_times_in;
    this->exact_ccp_rmesh_times = exact_ccp_rmesh_times_in;
    this->ewald_lambda = ewald_lambda_in;
    this->tail_mode = EwaldVqDetail::parse_tail_mode(GlobalC::exx_info.info_ri.ewald_tail_check);
    this->coulomb_param = coulomb_param_in;
    this->ewald_dimension = ewald_dimension_in;
    if (this->ewald_dimension != 2 && this->ewald_dimension != 3)
    {
        ModuleBase::WARNING_QUIT("Ewald_Vq::init", "ewald_dimension must be 2 or 3");
    }

    this->g_lcaos = this->init_gauss(lcaos_in);
    this->g_abfs = this->init_gauss(abfs_in);
    this->g_abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->g_abfs,
                                                        this->coulomb_param,
                                                        this->exact_ccp_rmesh_times);
    this->multipole = Exx_Abfs::Construct_Orbs::get_multipole(abfs_in);
    this->abfs = abfs_in;
    this->abfs_old_to_new = abfs_old_to_new;
    this->bare_multipole_scale = 0.0;
    const auto fock_params = this->coulomb_param.find(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
    if (fock_params != this->coulomb_param.end())
    {
        for (const auto& param: fock_params->second)
        {
            const auto alpha = param.find("alpha");
            if (alpha != param.end() && !alpha->second.empty())
            {
                this->bare_multipole_scale += std::stod(alpha->second);
            }
        }
    }
    this->lcaos_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(lcaos_in);
    this->abfs_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(abfs_in);
    this->g_lcaos_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(this->g_lcaos);
    this->g_abfs_ccp_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(this->g_abfs_ccp);

    const ModuleBase::Element_Basis_Index::Range range_abfs = ModuleBase::Element_Basis_Index::construct_range(abfs_in);
    this->index_abfs = ModuleBase::Element_Basis_Index::construct_index(range_abfs, abfs_old_to_new);

    this->cv
        .set_orbitals(ucell,
                      orb,
                      this->g_lcaos,
                      this->g_abfs,
                      this->g_abfs_ccp,
                      kmesh_times_in,
                      MGT_in,
                      false,
                      abfs_old_to_new);
    this->gaunt.create(MGT_in->Gaunt_Coefficients.getBound1(),
                       MGT_in->Gaunt_Coefficients.getBound2(),
                       MGT_in->Gaunt_Coefficients.getBound3());
    this->gaunt = MGT_in->Gaunt_Coefficients;

    this->atoms_vec.resize(ucell.nat);
    std::iota(this->atoms_vec.begin(), this->atoms_vec.end(), 0);
    this->nmp = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};
    if (this->ewald_dimension == 2)
    {
        const double eps = 1e-10;
        const ModuleBase::Vector3<double> a1 = ucell.a1;
        const ModuleBase::Vector3<double> a2 = ucell.a2;
        const ModuleBase::Vector3<double> a3 = ucell.a3;
        const double a3_norm = a3.norm();
        if (this->nmp[2] != 1)
        {
            ModuleBase::WARNING_QUIT("Ewald_Vq::init", "2D Ewald requires nkz = 1");
        }
        for (size_t ik = 0; ik != this->p_kv->kvec_c_full.size(); ++ik)
        {
            if (std::abs(this->p_kv->kvec_c_full[ik].z) > eps)
            {
                ModuleBase::WARNING_QUIT("Ewald_Vq::init", "2D Ewald requires zero k-point offset along z");
            }
        }
        if (std::abs(a1.z) > eps || std::abs(a2.z) > eps || std::abs(a3.x) > eps || std::abs(a3.y) > eps
            || std::abs(a1 * a3) > eps * a1.norm() * a3_norm
            || std::abs(a2 * a3) > eps * a2.norm() * a3_norm)
        {
            ModuleBase::WARNING_QUIT("Ewald_Vq::init", "2D Ewald requires a slab cell with the non-periodic direction along z");
        }
        if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
        {
            ModuleBase::WARNING_QUIT("Ewald_Vq::init", "2D Ewald force/stress derivatives are not implemented");
        }
    }

    ModuleBase::timer::end("Ewald_Vq", "init");
}

template <typename Tdata>
void Ewald_Vq<Tdata>::init_ions(const UnitCell& ucell, const std::array<Tcell, Ndim>& period_Vs_NAO)
{
    ModuleBase::TITLE("Ewald_Vq", "init_ions");
    ModuleBase::timer::start("Ewald_Vq", "init_ions");

    // The bare and Gaussian real-space maps are subtracted locally below.
    // Enumerate both with the bare-Coulomb period so every common key has the
    // same MPI owner; the Gaussian cutoff is still applied when building it.
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = EwaldVqDetail::distribute_common_realspace_tasks(this->mpi_comm,
                                                          this->atoms_vec,
                                                          period_Vs_NAO);

    this->list_A0 = list_As_Vs.first;
    this->list_A1 = list_As_Vs.second[0];

    const std::array<int, 1> Nks = {this->nks0};
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TK>>>> list_As_Vq
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, this->atoms_vec, Nks, 2, false);
    this->list_A0_k = list_As_Vq.first;
    this->list_A1_k = list_As_Vq.second[0];

    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, this->atoms_vec, period_Vs_NAO, 2, false);
    this->list_A0_pair_R = list_As_Vs_atoms.first;
    this->list_A1_pair_R = list_As_Vs_atoms.second[0];

    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms_period
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, this->atoms_vec, this->nmp, 2, false);
    this->list_A0_pair_R_period = list_As_Vs_atoms_period.first;
    this->list_A1_pair_R_period = list_As_Vs_atoms_period.second[0];

    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TK>>>> list_As_Vq_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, this->atoms_vec, Nks, 2, false);
    this->list_A0_pair_k = list_As_Vq_atoms.first;
    this->list_A1_pair_k = list_As_Vq_atoms.second[0];

    for (size_t ik = 0; ik != this->nks0; ++ik)
        this->kvec_c[ik] = this->p_kv->kvec_c_full[ik];

    std::vector<ModuleBase::Vector3<double>> neg_kvec(this->nks0);
    std::transform(this->kvec_c.begin(),
                   this->kvec_c.end(),
                   neg_kvec.begin(),
                   [](ModuleBase::Vector3<double>& vec) -> ModuleBase::Vector3<double> { return -vec; });
    this->gaussian_abfs.init(ucell, 2 * this->abfs_Lmax + 1, neg_kvec, ucell.G, this->ewald_lambda);

    ModuleBase::timer::end("Ewald_Vq", "init_ions");
}

template <typename Tdata>
double Ewald_Vq<Tdata>::get_singular_chi(const UnitCell& ucell, const std::vector<std::map<std::string,std::string>>& param_list, const double& qdiv)
{
    ModuleBase::TITLE("Ewald_Vq", "get_singular_chi");
    ModuleBase::timer::start("Ewald_Vq", "get_singular_chi");

    double chi = 0.0;
    for(const auto &param : param_list)
	{
        if(param.at("singularity_correction") == "carrier")
		{
            if (this->ewald_dimension == 2)
            {
                ModuleBase::WARNING_QUIT("Ewald_Vq::get_singular_chi", "2D Ewald supports massidda correction only");
            }
            chi = Singular_Value::cal_carrier(ucell, this->kvec_c, qdiv, 100, 30, 1e-6, 3);
        }
        else if(param.at("singularity_correction") == "massidda")
		{
            if (this->ewald_dimension == 2)
            {
                chi = Singular_Value::cal_massidda_2d(ucell, this->nmp, 1.0, 5, 1e-4);
            }
            else
            {
                chi = Singular_Value::cal_massidda(ucell, this->nmp, qdiv, 1, 5, 1e-4);
            }
        }
        else
        {            
            throw std::domain_error(std::string(__FILE__) + " line " + std::to_string(__LINE__)
                                    + ": singularity_correction must be carrier or massidda");
        }
    }

    ModuleBase::timer::end("Ewald_Vq", "get_singular_chi");
    return chi;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_gauss(const UnitCell& ucell, const std::vector<TA>& list_A0, const std::vector<TAC>& list_A1)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vs_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_Vs_gauss");

    std::map<std::string, bool> flags = {{"writable_Vws", true}};
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_gauss = this->cv.cal_Vs(ucell, list_A0, list_A1, flags);
    this->cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_gauss);

    ModuleBase::timer::end("Ewald_Vq", "cal_Vs_gauss");
    return Vs_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVs_gauss(const UnitCell& ucell, const std::vector<TA>& list_A0, const std::vector<TAC>& list_A1)
    -> std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVs_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVs_gauss");

    std::map<std::string, bool> flags = {{"writable_dVws", true}};

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_gauss
        = this->cv.cal_dVs(ucell, list_A0, list_A1, flags);
    this->cv.dVws = LRI_CV_Tools::get_dCVws(ucell, dVs_gauss);

    ModuleBase::timer::end("Ewald_Vq", "cal_dVs_gauss");
    return dVs_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_minus_gauss(const UnitCell& ucell,
                                         const std::vector<TA>& list_A0,
                                         const std::vector<TAC>& list_A1,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in,
                                         const bool truncate_analytic_tail)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vs_minus_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_Vs_minus_gauss");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_gauss = this->cal_Vs_gauss(ucell, list_A0, list_A1);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_minus_gauss;
    if (truncate_analytic_tail)
    {
        Vs_minus_gauss = this->set_Vs_minus_gauss_common_candidates(ucell,
                                                                    list_A0,
                                                                    list_A1,
                                                                    Vs_in,
                                                                    Vs_gauss);
    }
    else
    {
        Vs_minus_gauss = this->set_Vs_dVs_minus_gauss(ucell,
                                                       list_A0,
                                                       list_A1,
                                                       Vs_in,
                                                       Vs_gauss);
    }

    ModuleBase::timer::end("Ewald_Vq", "cal_Vs_minus_gauss");
    return Vs_minus_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVs_minus_gauss(const UnitCell& ucell,
                                          const std::vector<TA>& list_A0,
                                          const std::vector<TAC>& list_A1,
                                          std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_in)
    -> std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVs_minus_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVs_minus_gauss");

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_gauss
        = this->cal_dVs_gauss(ucell, list_A0, list_A1);
    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_minus_gauss
        = this->set_Vs_dVs_minus_gauss(ucell, list_A0, list_A1, dVs_in, dVs_gauss);

    ModuleBase::timer::end("Ewald_Vq", "cal_dVs_minus_gauss");
    return dVs_minus_gauss;
}

template <typename Tdata>
double Ewald_Vq<Tdata>::cal_V_Rcut(const int it0, const int it1)
{
    return this->g_abfs_ccp_rcut.at(it0) + this->g_lcaos_rcut.at(it1);
}

template <typename Tdata>
double Ewald_Vq<Tdata>::get_Rcut_max(const int it0, const int it1)
{
    double lcaos_rmax = this->lcaos_rcut.at(it0) * this->ccp_rmesh_times + this->lcaos_rcut.at(it1);
    double g_lcaos_rmax = this->g_lcaos_rcut.at(it0) * this->ccp_rmesh_times + this->g_lcaos_rcut.at(it1);
    return std::min(lcaos_rmax, g_lcaos_rmax);
}

template <typename Tdata>
template <typename Tresult>
auto Ewald_Vq<Tdata>::project_Vs_dVs_gauss(const UnitCell& ucell,
                                           const std::vector<TA>& list_A0,
                                           const std::vector<TAC>& list_A1,
                                           std::map<TA, std::map<TAC, Tresult>>& Vs_dVs_gauss_in)
    -> std::map<TA, std::map<TAC, Tresult>>
{
    using Tin_convert = typename LRI_CV_Tools::TinType<Tresult>::type;
    std::map<TA, std::map<TAC, Tresult>> pVs_dVs_gauss;
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (size_t i0 = 0; i0 < list_A0.size(); ++i0)
    {
        for (size_t i1 = 0; i1 < list_A1.size(); ++i1)
        {
            const TA iat0 = list_A0[i0];
            const int it0 = ucell.iat2it[iat0];
            const int ia0 = ucell.iat2ia[iat0];
            const TA iat1 = list_A1[i1].first;
            const int it1 = ucell.iat2it[iat1];
            const int ia1 = ucell.iat2ia[iat1];
            const TC& cell1 = list_A1[i1].second;

            const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];
            const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];

            const double Rcut = std::min(this->cal_V_Rcut(it0, it1), this->cal_V_Rcut(it1, it0));
            const Abfs::Vector3_Order<double> R_delta
                = -tau0 + tau1 + (RI_Util::array3_to_Vector3(cell1) * ucell.latvec);
            if (R_delta.norm() * ucell.lat0 < Rcut)
            {
                const size_t size0 = this->index_abfs[it0].count_size;
                const size_t size1 = this->index_abfs[it1].count_size;
                Tresult data;
                LRI_CV_Tools::init_elem(data, size0, size1);

                // pA * pB * V(R)_gauss
                for (int l0 = 0; l0 != this->g_abfs_ccp.at(it0).size(); ++l0)
                {
                    for (int l1 = 0; l1 != this->g_abfs.at(it1).size(); ++l1)
                    {
                        for (size_t n0 = 0; n0 != this->g_abfs_ccp.at(it0).at(l0).size(); ++n0)
                        {
                            const double pA = this->multipole.at(it0).at(l0).at(n0);
                            for (size_t n1 = 0; n1 != this->g_abfs.at(it1).at(l1).size(); ++n1)
                            {
                                const double pB = this->multipole.at(it1).at(l1).at(n1);
                                Tin_convert pp = RI::Global_Func::convert<Tin_convert>(pA * pB);
                                for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                                {
                                    for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                    {
                                        const size_t index0 = this->index_abfs.at(it0).at(l0).at(n0).at(m0);
                                        const size_t index1 = this->index_abfs.at(it1).at(l1).at(n1).at(m1);

                                        LRI_CV_Tools::add_elem(data,
                                                               index0,
                                                               index1,
                                                               Vs_dVs_gauss_in.at(list_A0[i0]).at(list_A1[i1]),
                                                               index0,
                                                               index1,
                                                               pp);
                                    }
                                }
                            }
                        }
                    }
                }
#pragma omp critical(Ewald_Vq_set_Vs_dVs_minus_gauss)
                pVs_dVs_gauss[list_A0[i0]][list_A1[i1]] = data;
            }
        }
    }

    return pVs_dVs_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_analytic_multipole_Vs(
    const UnitCell& ucell,
    const std::vector<EwaldVqDetail::TailKey>& keys)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> multipole_blocks;
    if (keys.empty()) return multipole_blocks;
    if (std::abs(this->bare_multipole_scale) <= std::numeric_limits<double>::epsilon())
    {
        throw std::invalid_argument("Failed to determine the Ewald Fock alpha for analytic multipole evaluation.");
    }

    Moment_abfs<Tdata> moment_abfs(GlobalC::exx_info.info_ri);
    moment_abfs.cal_multipole(this->abfs);
    const std::vector<double> zero_cutoff(this->abfs.size(), 0.0);
    const auto grouped = EwaldVqDetail::group_tail_keys_by_atom(keys);
    for (typename std::map<int, std::vector<TAC>>::const_iterator grouped_it = grouped.begin();
         grouped_it != grouped.end();
         ++grouped_it)
    {
        std::vector<TAC> nonzero_keys;
        const int iat0 = grouped_it->first;
        const int it0 = ucell.iat2it[iat0];
        const int ia0 = ucell.iat2ia[iat0];
        const auto tau0 = ucell.atoms[it0].tau[ia0];
        for (const TAC& jr: grouped_it->second)
        {
            const int it1 = ucell.iat2it[jr.first];
            const int ia1 = ucell.iat2ia[jr.first];
            const auto tau1 = ucell.atoms[it1].tau[ia1];
            const auto delta = -tau0 + tau1 + (RI_Util::array3_to_Vector3(jr.second) * ucell.latvec);
            if (delta.norm() * ucell.lat0 > 1e-12) nonzero_keys.push_back(jr);
        }
        if (nonzero_keys.empty()) continue;

        std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> part;
        moment_abfs.cal_VR(ucell,
                           this->abfs,
                           std::make_pair(std::vector<TA>{iat0}, std::vector<std::vector<TAC>>{nonzero_keys}),
                           zero_cutoff,
                           0.0,
                           this->cv,
                           part,
                           this->abfs_old_to_new,
                           false,
                           false,
                           false,
                           true,
                           this->bare_multipole_scale,
                           true);
        for (typename std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>::const_iterator part_it = part.begin();
             part_it != part.end();
             ++part_it)
        {
            multipole_blocks[part_it->first].insert(part_it->second.begin(), part_it->second.end());
        }
    }
    return multipole_blocks;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::set_Vs_minus_gauss_common_candidates(
    const UnitCell& ucell,
    const std::vector<TA>& list_A0,
    const std::vector<TAC>& list_A1,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_gauss_in)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    const auto projected_gaussian = this->project_Vs_dVs_gauss(ucell, list_A0, list_A1, Vs_gauss_in);
    const auto find_tensor = [](const auto& tensor_map, const EwaldVqDetail::TailKey& key)
        -> const RI::Tensor<Tdata>*
    {
        const auto outer = tensor_map.find(key.atom_i);
        if (outer == tensor_map.end()) return nullptr;
        const TAC jr{key.atom_j, key.cell};
        const auto inner = outer->second.find(jr);
        return inner == outer->second.end() ? nullptr : &inner->second;
    };

    std::vector<EwaldVqDetail::TailKey> multipole_keys;
    bool local_complete = true;
    EwaldVqDetail::TailKey first_incomplete_key{};
    double first_incomplete_distance = 0.0;
    double first_incomplete_bare_support = 0.0;
    double first_incomplete_gaussian_support = 0.0;
    bool first_incomplete_has_bare = false;
    bool first_incomplete_has_gaussian = false;
    for (const TA iat0: list_A0)
    {
        const int it0 = ucell.iat2it[iat0];
        const int ia0 = ucell.iat2ia[iat0];
        const auto tau0 = ucell.atoms[it0].tau[ia0];
        for (const TAC& jr: list_A1)
        {
            const EwaldVqDetail::TailKey key{iat0, jr.first, jr.second};
            const int it1 = ucell.iat2it[jr.first];
            const int ia1 = ucell.iat2ia[jr.first];
            const auto tau1 = ucell.atoms[it1].tau[ia1];
            const auto delta = -tau0 + tau1 + (RI_Util::array3_to_Vector3(jr.second) * ucell.latvec);
            const double distance = delta.norm() * ucell.lat0;
            const double bare_support = this->abfs_rcut[it0] + this->abfs_rcut[it1];
            const double gaussian_support
                = std::min(this->cal_V_Rcut(it0, it1), this->cal_V_Rcut(it1, it0));
            const auto plan = EwaldVqDetail::plan_tail_block(find_tensor(Vs_in, key) != nullptr,
                                                             find_tensor(projected_gaussian, key) != nullptr,
                                                             distance,
                                                             bare_support,
                                                             gaussian_support);
            if (plan.coverage == EwaldVqDetail::TailBlockCoverage::Incomplete)
            {
                if (local_complete)
                {
                    first_incomplete_key = key;
                    first_incomplete_distance = distance;
                    first_incomplete_bare_support = bare_support;
                    first_incomplete_gaussian_support = gaussian_support;
                    first_incomplete_has_bare = find_tensor(Vs_in, key) != nullptr;
                    first_incomplete_has_gaussian = find_tensor(projected_gaussian, key) != nullptr;
                }
                local_complete = false;
            }
            else if (plan.use_bare_multipole || plan.use_gaussian_multipole)
            {
                multipole_keys.push_back(key);
            }
        }
    }
    int complete = local_complete ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &complete, 1, MPI_INT, MPI_MIN, this->mpi_comm);
    if (complete == 0)
    {
        if (!local_complete)
        {
            int mpi_rank = 0;
            MPI_Comm_rank(this->mpi_comm, &mpi_rank);
            std::cerr << "Ewald candidate coverage failure on rank " << mpi_rank << ": I,J,R="
                      << first_incomplete_key.atom_i << ' ' << first_incomplete_key.atom_j << ' '
                      << first_incomplete_key.cell[0] << ' ' << first_incomplete_key.cell[1] << ' '
                      << first_incomplete_key.cell[2] << ", distance=" << first_incomplete_distance
                      << ", bare_support=" << first_incomplete_bare_support
                      << ", Gaussian_support=" << first_incomplete_gaussian_support
                      << ", has_bare=" << (first_incomplete_has_bare ? "true" : "false")
                      << ", has_Gaussian=" << (first_incomplete_has_gaussian ? "true" : "false")
                      << std::endl;
        }
        throw std::runtime_error("An Ewald candidate block is missing inside its compact support.");
    }

    const auto multipole_blocks = this->cal_analytic_multipole_Vs(ucell, multipole_keys);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> difference;
    for (const TA iat0: list_A0)
    {
        const int it0 = ucell.iat2it[iat0];
        const int ia0 = ucell.iat2ia[iat0];
        const auto tau0 = ucell.atoms[it0].tau[ia0];
        for (const TAC& jr: list_A1)
        {
            const EwaldVqDetail::TailKey key{iat0, jr.first, jr.second};
            const int it1 = ucell.iat2it[jr.first];
            const int ia1 = ucell.iat2ia[jr.first];
            const auto tau1 = ucell.atoms[it1].tau[ia1];
            const auto delta = -tau0 + tau1 + (RI_Util::array3_to_Vector3(jr.second) * ucell.latvec);
            const double distance = delta.norm() * ucell.lat0;
            const double bare_support = this->abfs_rcut[it0] + this->abfs_rcut[it1];
            const double gaussian_support
                = std::min(this->cal_V_Rcut(it0, it1), this->cal_V_Rcut(it1, it0));
            const RI::Tensor<Tdata>* bare = find_tensor(Vs_in, key);
            const RI::Tensor<Tdata>* gaussian = find_tensor(projected_gaussian, key);
            const auto plan = EwaldVqDetail::plan_tail_block(bare != nullptr,
                                                             gaussian != nullptr,
                                                             distance,
                                                             bare_support,
                                                             gaussian_support);
            if (plan.coverage == EwaldVqDetail::TailBlockCoverage::AnalyticCancellation) continue;
            const RI::Tensor<Tdata>* analytic = find_tensor(multipole_blocks, key);
            if (plan.use_bare_multipole) bare = analytic;
            if (plan.use_gaussian_multipole) gaussian = analytic;
            if (bare == nullptr || gaussian == nullptr)
            {
                throw std::runtime_error("Failed to evaluate an Ewald analytic multipole block.");
            }
            difference[iat0][jr] = *bare - *gaussian;
        }
    }
    return difference;
}

template <typename Tdata>
template <typename Tresult>
auto Ewald_Vq<Tdata>::set_Vs_dVs_minus_gauss(const UnitCell& ucell,
                                             const std::vector<TA>& list_A0,
                                             const std::vector<TAC>& list_A1,
                                             std::map<TA, std::map<TAC, Tresult>>& Vs_dVs_in,
                                             std::map<TA, std::map<TAC, Tresult>>& Vs_dVs_gauss_in)
    -> std::map<TA, std::map<TAC, Tresult>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vs_dVs_minus_gauss");
    ModuleBase::timer::tick("Ewald_Vq", "set_Vs_dVs_minus_gauss");

    auto pVs_dVs_gauss = this->project_Vs_dVs_gauss(ucell, list_A0, list_A1, Vs_dVs_gauss_in);
    auto difference = LRI_CV_Tools::minus_common_keys(Vs_dVs_in, pVs_dVs_gauss);
    ModuleBase::timer::tick("Ewald_Vq", "set_Vs_dVs_minus_gauss");
    return difference;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vq_gauss(const UnitCell& ucell,
                                   const std::vector<TA>& list_A0_k,
                                   const std::vector<TAK>& list_A1_k,
                                   const double& chi,
                                   const int& shift_for_mpi)
    -> std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vq_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_Vq_gauss");

    T_func_DPget_Vq_dVq<RI::Tensor<std::complex<double>>> func_DPget_Vq;
    if (this->ewald_dimension == 2)
    {
        func_DPget_Vq = std::bind(&Gaussian_Abfs::get_Vq_2d,
                                  &this->gaussian_abfs,
                                  std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3,
                                  chi,
                                  std::placeholders::_4,
                                  this->gaunt);
    }
    else
    {
        func_DPget_Vq = std::bind(&Gaussian_Abfs::get_Vq,
                                  &this->gaussian_abfs,
                                  std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3,
                                  chi,
                                  std::placeholders::_4,
                                  this->gaunt);
    }
    auto Vq_gauss = this->set_Vq_dVq_gauss(ucell, list_A0_k, list_A1_k, shift_for_mpi, func_DPget_Vq);

    ModuleBase::timer::end("Ewald_Vq", "cal_Vq_gauss");
    return Vq_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVq_gauss(const UnitCell& ucell,
                                    const std::vector<TA>& list_A0_k,
                                    const std::vector<TAK>& list_A1_k,
                                    const double& chi,
                                    const int& shift_for_mpi)
    -> std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVq_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVq_gauss");

    if (this->ewald_dimension == 2)
    {
        ModuleBase::WARNING_QUIT("Ewald_Vq::cal_dVq_gauss", "2D Ewald force/stress derivatives are not implemented");
    }

    const T_func_DPget_Vq_dVq<std::array<RI::Tensor<std::complex<double>>, Ndim>> func_DPget_dVq
        = std::bind(&Gaussian_Abfs::get_dVq,
                    &this->gaussian_abfs,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3,
                    chi,
                    std::placeholders::_4,
                    this->gaunt);

    using namespace RI::Array_Operator;
    std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>> dVq_gauss;
    auto res = this->set_Vq_dVq_gauss(ucell, list_A0_k, list_A1_k, shift_for_mpi, func_DPget_dVq);

    for (size_t i0 = 0; i0 < list_A0_k.size(); ++i0)
    {
        const TA iat0 = list_A0_k[i0];
        const int it0 = ucell.iat2it[iat0];
        for (size_t i1 = 0; i1 < list_A1_k.size(); ++i1)
        {
            const TA iat1 = list_A1_k[i1].first;
            const int it1 = ucell.iat2it[iat1];
            if (iat0 != iat1)
            {
                const int ik = list_A1_k[i1].second[0];
                const TAK index0 = std::make_pair(iat1, TK{ik});
                dVq_gauss[iat0][index0] = -res.at(list_A0_k[i0]).at(list_A1_k[i1]);
                const TAK index1 = std::make_pair(iat0, TK{ik});
                dVq_gauss[iat1][index1] = res.at(list_A0_k[i0]).at(list_A1_k[i1]);
            }
            else
                dVq_gauss[list_A0_k[i0]][list_A1_k[i1]] = res.at(list_A0_k[i0]).at(list_A1_k[i1]);
        }
    }

    ModuleBase::timer::end("Ewald_Vq", "cal_dVq_gauss");
    return dVq_gauss;
}

template <typename Tdata>
template <typename Tresult>
auto Ewald_Vq<Tdata>::set_Vq_dVq_gauss(const UnitCell& ucell,
                                       const std::vector<TA>& list_A0_k,
                                       const std::vector<TAK>& list_A1_k,
                                       const int& shift_for_mpi,
                                       const T_func_DPget_Vq_dVq<Tresult>& func_DPget_Vq_dVq)
    -> std::map<TA, std::map<TAK, Tresult>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vq_dVq_gauss");
    ModuleBase::timer::start("Ewald_Vq", "set_Vq_dVq_gauss");

    std::map<TA, std::map<TAK, Tresult>> Vq_dVq_gauss_out;
    for(const auto &param_list : this->coulomb_param)
	{
        std::complex<double> alpha;
        for(const auto &param : param_list.second)
        {
            alpha = std::complex<double>(std::stod(param.at("alpha")), 0);
        }
#pragma omp parallel for collapse(2) schedule(dynamic)
        for (size_t i0 = 0; i0 < list_A0_k.size(); ++i0)
        {
            for (size_t i1 = 0; i1 < list_A1_k.size(); ++i1)
            {
                const TA iat0 = list_A0_k[i0];
                const int it0 = ucell.iat2it[iat0];
                const int ia0 = ucell.iat2ia[iat0];
                const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];

                const TA iat1 = list_A1_k[i1].first;
                const int it1 = ucell.iat2it[iat1];
                const int ia1 = ucell.iat2ia[iat1];
                const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
                const size_t ik = list_A1_k[i1].second[0] + shift_for_mpi;

                const ModuleBase::Vector3<double> tau = tau0 - tau1;
                auto data
                    = func_DPget_Vq_dVq(this->g_abfs_ccp.at(it0).size() - 1, this->g_abfs.at(it1).size() - 1, ik, tau);

#pragma omp critical(Ewald_Vq_set_Vq_dVq_gauss)
                Vq_dVq_gauss_out[list_A0_k[i0]][list_A1_k[i1]] = LRI_CV_Tools::mul2(alpha, data);
            }
        }
    }

    ModuleBase::timer::end("Ewald_Vq", "set_Vq_dVq_gauss");
    return Vq_dVq_gauss_out;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vq_minus_gauss(const UnitCell& ucell,
                                         const std::vector<TA>& list_A0,
                                         const std::vector<TAC>& list_A1,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_minus_gauss)
    -> std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vq_minus_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_Vq_minus_gauss");

    auto Vq_minus_gauss
        = this->set_Vq_dVq_minus_gauss<RI::Tensor<std::complex<double>>>(ucell, list_A0, list_A1, Vs_minus_gauss);

    ModuleBase::timer::end("Ewald_Vq", "cal_Vq_minus_gauss");
    return Vq_minus_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVq_minus_gauss(
    const UnitCell& ucell,
    const std::vector<TA>& list_A0,
    const std::vector<TAC>& list_A1,
    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_minus_gauss)
    -> std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVq_minus_gauss");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVq_minus_gauss");

    auto dVq_minus_gauss
        = this->set_Vq_dVq_minus_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>>(ucell,
                                                                                           list_A0,
                                                                                           list_A1,
                                                                                           dVs_minus_gauss);

    ModuleBase::timer::end("Ewald_Vq", "cal_dVq_minus_gauss");
    return dVq_minus_gauss;
}

template <typename Tdata>
template <typename Tout, typename Tin>
auto Ewald_Vq<Tdata>::set_Vq_dVq_minus_gauss(const UnitCell& ucell,
                                             const std::vector<TA>& list_A0,
                                             const std::vector<TAC>& list_A1,
                                             std::map<TA, std::map<TAC, Tin>>& Vs_dVs_minus_gauss)
    -> std::map<TA, std::map<TAK, Tout>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vq_dVq_minus_gauss");
    ModuleBase::timer::start("Ewald_Vq", "set_Vq_dVq_minus_gauss");

    using namespace RI::Array_Operator;
    using Tin_convert = typename LRI_CV_Tools::TinType<Tout>::type;
    std::map<TA, std::map<TAK, Tout>> datas;

    // auto start = std::chrono::system_clock::now();

#pragma omp parallel
    {
        std::map<TA, std::map<TAK, Tout>> local_datas;

#pragma omp for schedule(dynamic) nowait
        for (size_t ik = 0; ik != this->nks0; ++ik)
        {
            for (size_t i0 = 0; i0 < list_A0.size(); ++i0)
            {
                for (size_t i1 = 0; i1 < list_A1.size(); ++i1)
                {
                    const TA iat0 = list_A0[i0];
                    const int it0 = ucell.iat2it[iat0];
                    const int ia0 = ucell.iat2ia[iat0];
                    const TA iat1 = list_A1[i1].first;
                    const int it1 = ucell.iat2it[iat1];
                    const int ia1 = ucell.iat2ia[iat1];
                    const TC& cell1 = list_A1[i1].second;

                    const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];
                    const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
                    const ModuleBase::Vector3<double> R_delta
                        = -tau0 + tau1 + (RI_Util::array3_to_Vector3(cell1) * ucell.latvec);
                    const auto atom_blocks = Vs_dVs_minus_gauss.find(iat0);
                    if (atom_blocks == Vs_dVs_minus_gauss.end()) continue;
                    const auto block = atom_blocks->second.find(list_A1[i1]);
                    if (block == atom_blocks->second.end()) continue;

                    const double Rcut = std::min(this->get_Rcut_max(it0, it1),
                                                 this->get_Rcut_max(it1, it0));
                    if (EwaldVqDetail::uses_adaptive_realspace_range(this->tail_mode)
                        || R_delta.norm() * ucell.lat0 < Rcut)
                    {
                        const std::complex<double> phase
                            = std::exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT
                                       * (this->kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1) * ucell.latvec)));

                        Tin block_value = block->second;
                        Tout Vs_dVs_tmp = LRI_CV_Tools::mul2(
                            phase,
                            LRI_CV_Tools::convert<Tin_convert>(std::move(block_value)));

                        const TAK index = std::make_pair(iat1, TK{static_cast<int>(ik)});
                        if (!LRI_CV_Tools::exist(local_datas[iat0][index]))
                            local_datas[iat0][index] = Vs_dVs_tmp;
                        else
                            local_datas[iat0][index] = local_datas.at(iat0).at(index) + Vs_dVs_tmp;
                    }
                }
            }
        }

#pragma omp critical(Ewald_Vq_set_Vq_dVq_minus_gauss)
        {
            for (auto it0 = local_datas.begin(); it0 != local_datas.end(); ++it0)
            {
                const TA& key0 = it0->first;
                std::map<TAK, Tout>& map1 = it0->second;
                for (auto it1 = map1.begin(); it1 != map1.end(); ++it1)
                {
                    const TAK& key1 = it1->first;
                    Tout& value = it1->second;

                    if (!LRI_CV_Tools::exist(datas[key0][key1]))
                        datas[key0][key1] = value;
                    else
                        datas[key0][key1] = datas.at(key0).at(key1) + value;
                }
            }
        }
    }

    // auto end = std::chrono::system_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end
    // - start); std::cout << "set_Vq_dVq_minus_gauss Time: "
    //           << double(duration.count()) *
    //           std::chrono::microseconds::period::num
    //                  / std::chrono::microseconds::period::den
    //           << " s" << std::endl;
    ModuleBase::timer::end("Ewald_Vq", "set_Vq_dVq_minus_gauss");
    return datas;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vq(const UnitCell& ucell,
                             const double& chi,
                             std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in)
    -> std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vq");
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_minus_gauss = this->cal_Vs_minus_gauss(ucell,
                                                                                             this->list_A0,
                                                                                             this->list_A1,
                                                                                             Vs_in,
                                                                                             this->tail_mode == EwaldVqDetail::TailMode::Enlarge
                                                                                                 || this->tail_mode == EwaldVqDetail::TailMode::Strict); //{ia0, {ia1, R}}
    const T_func_DPcal_Vq_dVq_minus_gauss<RI::Tensor<std::complex<double>>, RI::Tensor<Tdata>> func_cal_Vq_minus_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_Vq_minus_gauss,
                    this,
                    std::ref(ucell),
                    this->list_A0_pair_R,
                    this->list_A1_pair_R,
                    std::placeholders::_1);
    const T_func_DPcal_Vq_dVq_gauss<RI::Tensor<std::complex<double>>> func_cal_Vq_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_Vq_gauss,
                    this,
                    std::ref(ucell),
                    this->list_A0_k,
                    this->list_A1_k,
                    chi,
                    std::placeholders::_1);

    auto Vq = this->set_Vq_dVq(ucell,
                               this->list_A0_pair_k,
                               this->list_A1_pair_k,
	                               Vs_minus_gauss,
	                               func_cal_Vq_minus_gauss,
	                               func_cal_Vq_gauss); //{ia0, ia1}

    ModuleBase::timer::end("Ewald_Vq", "cal_Vq");
    return Vq;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVq(const UnitCell& ucell,
                              const double& chi,
                              std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_in)
    -> std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVq");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVq");

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_minus_gauss
        = this->cal_dVs_minus_gauss(ucell,
                                    this->list_A0,
                                    this->list_A1,
                                    dVs_in); //{ia0, {ia1, R}}
    const T_func_DPcal_Vq_dVq_minus_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>,
                                          std::array<RI::Tensor<Tdata>, Ndim>>
        func_cal_dVq_minus_gauss = std::bind(&Ewald_Vq<Tdata>::cal_dVq_minus_gauss,
                                             this,
                                             std::ref(ucell),
                                             this->list_A0_pair_R,
                                             this->list_A1_pair_R,
                                             std::placeholders::_1);
    const T_func_DPcal_Vq_dVq_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>> func_cal_dVq_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_dVq_gauss,
                    this,
                    std::ref(ucell),
                    this->list_A0_k,
                    this->list_A1_k,
                    chi,
                    std::placeholders::_1);

    auto dVq = this->set_Vq_dVq(ucell,
                                this->list_A0_pair_k,
                                this->list_A1_pair_k,
                                dVs_minus_gauss,
                                func_cal_dVq_minus_gauss,
                                func_cal_dVq_gauss); //{ia0, ia1}

    ModuleBase::timer::end("Ewald_Vq", "cal_dVq");
    return dVq;
}

template <typename Tdata>
template <typename Tout, typename Tin>
auto Ewald_Vq<Tdata>::set_Vq_dVq(const UnitCell& ucell,
                                 const std::vector<TA>& list_A0_pair_k,
                                 const std::vector<TAK>& list_A1_pair_k,
                                 std::map<TA, std::map<TAC, Tin>>& Vs_dVs_minus_gauss_in,
                                 const T_func_DPcal_Vq_dVq_minus_gauss<Tout, Tin>& func_cal_Vq_dVq_minus_gauss,
                                 const T_func_DPcal_Vq_dVq_gauss<Tout>& func_cal_Vq_dVq_gauss)
    -> std::map<TA, std::map<TAK, Tout>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vq_dVq");
    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq");

    using namespace RI::Array_Operator;
    using Tin_convert = typename LRI_CV_Tools::TinType<Tout>::type;
    std::map<TA, std::map<TAK, Tout>> Vq_dVq;
    const int shift_for_mpi = std::floor(this->nks0 / 2.0);

    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: this->list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: this->list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }

    std::map<TA, std::map<TAC, Tin>> Vs_dVs_minus_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_dVs_minus_gauss_in, atoms00, atoms01);
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_minus_gauss_local
        = func_cal_Vq_dVq_minus_gauss(Vs_dVs_minus_gauss); // partial sum over local R blocks

    // MPI: combine partial Fourier sums from different local R subsets onto the
    // target (I, J, k) distribution.
    std::set<TA> atoms10;
    std::set<TA> atoms11;
    for (const auto& I: this->list_A0_pair_k)
    {
        atoms10.insert(I);
    }
    for (const auto& JR: this->list_A1_pair_k)
    {
        atoms11.insert(JR.first);
    }
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_minus_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_dVq_minus_gauss_local, atoms10, atoms11);

    // MPI: {ia0, {ia1, k}} to {ia0, ia1}
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_gauss_out = func_cal_Vq_dVq_gauss(shift_for_mpi); //{ia0, {ia1, k}}
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_dVq_gauss_out, atoms10, atoms11); //{ia0, ia1}

#pragma omp parallel for collapse(2) schedule(dynamic)
    for (size_t i0 = 0; i0 < list_A0_pair_k.size(); ++i0)
    {
        for (size_t i1 = 0; i1 < list_A1_pair_k.size(); ++i1)
        {
            const TA iat0 = list_A0_pair_k[i0];
            const int it0 = ucell.iat2it[iat0];
            const TA iat1 = list_A1_pair_k[i1].first;
            const int it1 = ucell.iat2it[iat1];
            const int ik = list_A1_pair_k[i1].second[0] + shift_for_mpi;
            const TAK re_index = std::make_pair(iat1, std::array<int, 1>{ik});

            // check the Fourier transformed V(q)
            // whether ccp_rmesh_times * Rcut >= rIJ
            // skip some IJ pairs
            auto it_outer = Vq_dVq_minus_gauss.find(list_A0_pair_k[i0]);
            if (it_outer == Vq_dVq_minus_gauss.end())
                continue;

            auto it_inner = it_outer->second.find(re_index);
            if (it_inner == it_outer->second.end())
                continue;

            const size_t size0 = this->index_abfs[it0].count_size;
            const size_t size1 = this->index_abfs[it1].count_size;
            Tout data;
            LRI_CV_Tools::init_elem(data, size0, size1);
            for (int l0 = 0; l0 != this->g_abfs_ccp.at(it0).size(); ++l0)
            {
                for (int l1 = 0; l1 != this->g_abfs.at(it1).size(); ++l1)
                {
                    for (size_t n0 = 0; n0 != this->g_abfs_ccp.at(it0).at(l0).size(); ++n0)
                    {
                        const double pA = this->multipole.at(it0).at(l0).at(n0);
                        for (size_t n1 = 0; n1 != this->g_abfs.at(it1).at(l1).size(); ++n1)
                        {
                            const double pB = this->multipole.at(it1).at(l1).at(n1);
                            Tin_convert frac = RI::Global_Func::convert<Tin_convert>(pA * pB);
                            for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                            {
                                const size_t index0 = this->index_abfs.at(it0).at(l0).at(n0).at(m0);
                                const size_t lm0 = l0 * l0 + m0;
                                for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                {
                                    const size_t index1 = this->index_abfs.at(it1).at(l1).at(n1).at(m1);
                                    const size_t lm1 = l1 * l1 + m1;

                                    LRI_CV_Tools::add_elem(data,
                                                           index0,
                                                           index1,
                                                           Vq_dVq_gauss.at(list_A0_pair_k[i0]).at(list_A1_pair_k[i1]),
                                                           lm0,
                                                           lm1,
                                                           frac);
                                }
                            }
                        }
                    }
                }
            }

#pragma omp critical(Ewald_Vq_set_Vq_dVq)
            if (LRI_CV_Tools::exist(Vq_dVq_minus_gauss.at(list_A0_pair_k[i0]).at(re_index)))
                Vq_dVq[list_A0_pair_k[i0]][re_index] = Vq_dVq_minus_gauss.at(list_A0_pair_k[i0]).at(re_index) + data;
        }
    }

    ModuleBase::timer::end("Ewald_Vq", "set_Vq_dVq");
    return Vq_dVq;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs(const UnitCell& ucell,
                             const double& chi,
                             std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in) //{ia0, {ia1, R}}
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
	ModuleBase::TITLE("Ewald_Vq", "cal_Vs");
	ModuleBase::timer::tick("Ewald_Vq", "cal_Vs");

    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq = this->cal_Vq(ucell, chi, Vs_in);
    auto Vs = this->set_Vs_dVs<RI::Tensor<Tdata>>(ucell,
                                                  this->list_A0_pair_R_period,
                                                  this->list_A1_pair_R_period,
                                                  Vq); //{ia0, ia1}

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs");
    return Vs;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_split(const UnitCell& ucell,
                                   const double& chi,
                                   std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in)
    -> std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>,
                 std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vs_split");
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_split");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_minus_gauss = this->cal_Vs_minus_gauss(ucell,
                                                                                             this->list_A0,
                                                                                             this->list_A1,
                                                                                             Vs_in,
                                                                                             this->tail_mode == EwaldVqDetail::TailMode::Enlarge
                                                                                                 || this->tail_mode == EwaldVqDetail::TailMode::Strict);
    const int shift_for_mpi = std::floor(this->nks0 / 2.0);
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: this->list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: this->list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_minus_gauss_comm
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_minus_gauss, atoms00, atoms01);
    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq_short_local
        = this->cal_Vq_minus_gauss(ucell,
                                   this->list_A0_pair_R,
                                   this->list_A1_pair_R,
                                   Vs_minus_gauss_comm);

    std::set<TA> atoms10;
    std::set<TA> atoms11;
    for (const auto& I: this->list_A0_pair_k)
    {
        atoms10.insert(I);
    }
    for (const auto& JR: this->list_A1_pair_k)
    {
        atoms11.insert(JR.first);
    }
    auto Vq_short = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_short_local, atoms10, atoms11);
    auto Vq_long_out = this->cal_Vq_gauss(ucell,
                                          this->list_A0_k,
                                          this->list_A1_k,
                                          chi,
                                          shift_for_mpi);
    auto Vq_long_gauss = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_long_out, atoms10, atoms11);
    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq_long;
    using Tin_convert = typename LRI_CV_Tools::TinType<RI::Tensor<std::complex<double>>>::type;

#pragma omp parallel for collapse(2) schedule(dynamic)
    for (size_t i0 = 0; i0 < this->list_A0_pair_k.size(); ++i0)
    {
        for (size_t i1 = 0; i1 < this->list_A1_pair_k.size(); ++i1)
        {
            const TA iat0 = this->list_A0_pair_k[i0];
            const int it0 = ucell.iat2it[iat0];
            const TA iat1 = this->list_A1_pair_k[i1].first;
            const int it1 = ucell.iat2it[iat1];
            const int ik = this->list_A1_pair_k[i1].second[0] + shift_for_mpi;
            const TAK re_index = std::make_pair(iat1, std::array<int, 1>{ik});

            auto it_short_outer = Vq_short.find(iat0);
            if (it_short_outer == Vq_short.end())
            {
                continue;
            }
            if (it_short_outer->second.find(re_index) == it_short_outer->second.end())
            {
                continue;
            }
            auto it_long_outer = Vq_long_gauss.find(iat0);
            if (it_long_outer == Vq_long_gauss.end())
            {
                continue;
            }
            auto it_long_inner = it_long_outer->second.find(this->list_A1_pair_k[i1]);
            if (it_long_inner == it_long_outer->second.end())
            {
                continue;
            }

            RI::Tensor<std::complex<double>> data;
            LRI_CV_Tools::init_elem(data, this->index_abfs[it0].count_size, this->index_abfs[it1].count_size);
            for (int l0 = 0; l0 != this->g_abfs_ccp.at(it0).size(); ++l0)
            {
                for (int l1 = 0; l1 != this->g_abfs.at(it1).size(); ++l1)
                {
                    for (size_t n0 = 0; n0 != this->g_abfs_ccp.at(it0).at(l0).size(); ++n0)
                    {
                        const double pA = this->multipole.at(it0).at(l0).at(n0);
                        for (size_t n1 = 0; n1 != this->g_abfs.at(it1).at(l1).size(); ++n1)
                        {
                            const double pB = this->multipole.at(it1).at(l1).at(n1);
                            Tin_convert frac = RI::Global_Func::convert<Tin_convert>(pA * pB);
                            for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                            {
                                const size_t index0 = this->index_abfs.at(it0).at(l0).at(n0).at(m0);
                                const size_t lm0 = l0 * l0 + m0;
                                for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                {
                                    const size_t index1 = this->index_abfs.at(it1).at(l1).at(n1).at(m1);
                                    const size_t lm1 = l1 * l1 + m1;
                                    LRI_CV_Tools::add_elem(data,
                                                           index0,
                                                           index1,
                                                           it_long_inner->second,
                                                           lm0,
                                                           lm1,
                                                           frac);
                                }
                            }
                        }
                    }
                }
            }

#pragma omp critical(Ewald_Vq_cal_Vs_split)
            if (LRI_CV_Tools::exist(data))
            {
                Vq_long[iat0][re_index] = std::move(data);
            }
        }
    }

    auto Vs_short = this->set_Vs_dVs<RI::Tensor<Tdata>>(ucell,
                                                        this->list_A0_pair_R_period,
                                                        this->list_A1_pair_R_period,
                                                        Vq_short);
    auto Vs_long = this->set_Vs_dVs<RI::Tensor<Tdata>>(ucell,
                                                       this->list_A0_pair_R_period,
                                                       this->list_A1_pair_R_period,
                                                       Vq_long);

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_split");
    return {Vs_short, Vs_long};
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_serial_full(const UnitCell& ucell,
                                         const double& chi,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in_full,
                                         const std::array<Tcell, Ndim>& period_Vs_NAO)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    const MPI_Comm mpi_comm_saved = this->mpi_comm;
    const auto list_A0_saved = this->list_A0;
    const auto list_A1_saved = this->list_A1;
    const auto list_A0_k_saved = this->list_A0_k;
    const auto list_A1_k_saved = this->list_A1_k;
    const auto list_A0_pair_R_saved = this->list_A0_pair_R;
    const auto list_A1_pair_R_saved = this->list_A1_pair_R;
    const auto list_A0_pair_R_period_saved = this->list_A0_pair_R_period;
    const auto list_A1_pair_R_period_saved = this->list_A1_pair_R_period;
    const auto list_A0_pair_k_saved = this->list_A0_pair_k;
    const auto list_A1_pair_k_saved = this->list_A1_pair_k;

    const std::array<Tcell, Ndim>& period_Vs = period_Vs_NAO;
    const std::array<int, 1> Nks = {this->nks0};

    this->mpi_comm = MPI_COMM_SELF;
    this->list_A0 = this->atoms_vec;
    this->list_A1 = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs);
    this->list_A0_k = this->atoms_vec;
    this->list_A1_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);
    this->list_A0_pair_R = this->atoms_vec;
    this->list_A1_pair_R = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs_NAO);
    this->list_A0_pair_R_period = this->atoms_vec;
    this->list_A1_pair_R_period = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, this->nmp);
    this->list_A0_pair_k = this->atoms_vec;
    this->list_A1_pair_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);

    auto Vs_full = this->cal_Vs(ucell, chi, Vs_in_full);

    this->mpi_comm = mpi_comm_saved;
    this->list_A0 = list_A0_saved;
    this->list_A1 = list_A1_saved;
    this->list_A0_k = list_A0_k_saved;
    this->list_A1_k = list_A1_k_saved;
    this->list_A0_pair_R = list_A0_pair_R_saved;
    this->list_A1_pair_R = list_A1_pair_R_saved;
    this->list_A0_pair_R_period = list_A0_pair_R_period_saved;
    this->list_A1_pair_R_period = list_A1_pair_R_period_saved;
    this->list_A0_pair_k = list_A0_pair_k_saved;
    this->list_A1_pair_k = list_A1_pair_k_saved;

    return Vs_full;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVs(
    const UnitCell& ucell,
    const double& chi,
    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_in) //{ia0, {ia1, R}}
    -> std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVs");
    ModuleBase::timer::start("Ewald_Vq", "cal_dVs");

    std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>> dVq
        = this->cal_dVq(ucell, chi, dVs_in);
    auto dVs = this->set_Vs_dVs<std::array<RI::Tensor<Tdata>, Ndim>>(ucell,
                                                                     this->list_A0_pair_R_period,
                                                                     this->list_A1_pair_R_period,
                                                                     dVq); //{ia0, ia1}

    ModuleBase::timer::end("Ewald_Vq", "cal_dVs");
    return dVs;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_realspace_difference(
    const UnitCell& ucell,
    const std::vector<TA>& list_A0,
    const std::vector<TAC>& list_A1,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in,
    const bool truncate_analytic_tail)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    return this->cal_Vs_minus_gauss(ucell,
                                    list_A0,
                                    list_A1,
                                    Vs_in,
                                    truncate_analytic_tail);
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_short_range_tail_stats(
    const UnitCell& ucell,
    const std::vector<EwaldVqDetail::TailKey>& local_keys,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_bare,
    const std::size_t local_production_blocks,
    const double local_reference_norm) -> EwaldVqDetail::TailStats
{
    using TailKey = EwaldVqDetail::TailKey;
    using TailBlockCoverage = EwaldVqDetail::TailBlockCoverage;

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> projected_gaussian;
    const auto grouped = EwaldVqDetail::group_tail_keys_by_atom(local_keys);
    for (std::map<int, std::vector<std::pair<int, EwaldVqDetail::TailCell>>>::const_iterator grouped_it
             = grouped.begin();
         grouped_it != grouped.end();
         ++grouped_it)
    {
        const int iat0 = grouped_it->first;
        const std::vector<std::pair<int, EwaldVqDetail::TailCell>>& list_A1 = grouped_it->second;
        const std::vector<TA> list_A0{iat0};
        auto gaussian = this->cal_Vs_gauss(ucell, list_A0, list_A1);
        auto projected = this->project_Vs_dVs_gauss(ucell, list_A0, list_A1, gaussian);
        for (typename std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>::const_iterator projected_it
                 = projected.begin();
             projected_it != projected.end();
             ++projected_it)
        {
            projected_gaussian[projected_it->first].insert(projected_it->second.begin(),
                                                            projected_it->second.end());
        }
    }

    const auto find_tensor = [](const auto& tensor_map, const TailKey& key)
        -> const RI::Tensor<Tdata>*
    {
        const auto outer = tensor_map.find(key.atom_i);
        if (outer == tensor_map.end()) return nullptr;
        const TAC jr{key.atom_j, key.cell};
        const auto inner = outer->second.find(jr);
        return inner == outer->second.end() ? nullptr : &inner->second;
    };

    std::map<TailKey, EwaldVqDetail::TailBlockPlan> plans;
    std::vector<TailKey> multipole_keys;
    for (const TailKey& key: local_keys)
    {
        const int it0 = ucell.iat2it[key.atom_i];
        const int ia0 = ucell.iat2ia[key.atom_i];
        const int it1 = ucell.iat2it[key.atom_j];
        const int ia1 = ucell.iat2ia[key.atom_j];
        const auto tau0 = ucell.atoms[it0].tau[ia0];
        const auto tau1 = ucell.atoms[it1].tau[ia1];
        const auto delta = -tau0 + tau1 + (RI_Util::array3_to_Vector3(key.cell) * ucell.latvec);
        const double distance = delta.norm() * ucell.lat0;
        const double bare_support = this->abfs_rcut[it0] + this->abfs_rcut[it1];
        const double gaussian_support
            = std::min(this->cal_V_Rcut(it0, it1), this->cal_V_Rcut(it1, it0));

        const RI::Tensor<Tdata>* bare = find_tensor(Vs_bare, key);
        const RI::Tensor<Tdata>* gaussian = find_tensor(projected_gaussian, key);
        const auto plan = EwaldVqDetail::plan_tail_block(bare != nullptr,
                                                         gaussian != nullptr,
                                                         distance,
                                                         bare_support,
                                                         gaussian_support);
        plans[key] = plan;
        if (plan.use_bare_multipole || plan.use_gaussian_multipole)
        {
            multipole_keys.push_back(key);
        }
    }
    const auto multipole_blocks = this->cal_analytic_multipole_Vs(ucell, multipole_keys);

    EwaldVqDetail::TailStats local_stats;
    local_stats.production_blocks = local_production_blocks;
    local_stats.reference_norm = local_reference_norm;
    std::map<TailKey, TailBlockCoverage> coverage;
    std::map<TailKey, RI::Tensor<Tdata>> differences;
    for (const TailKey& key: local_keys)
    {
        const int it0 = ucell.iat2it[key.atom_i];
        const int ia0 = ucell.iat2ia[key.atom_i];
        const int it1 = ucell.iat2it[key.atom_j];
        const int ia1 = ucell.iat2ia[key.atom_j];
        const auto tau0 = ucell.atoms[it0].tau[ia0];
        const auto tau1 = ucell.atoms[it1].tau[ia1];
        const auto delta = -tau0 + tau1 + (RI_Util::array3_to_Vector3(key.cell) * ucell.latvec);
        const double distance = delta.norm() * ucell.lat0;
        const auto plan = plans.at(key);
        coverage[key] = plan.coverage;

        const RI::Tensor<Tdata>* bare = find_tensor(Vs_bare, key);
        const RI::Tensor<Tdata>* gaussian = find_tensor(projected_gaussian, key);
        const RI::Tensor<Tdata>* analytic = find_tensor(multipole_blocks, key);
        if (plan.use_bare_multipole) bare = analytic;
        if (plan.use_gaussian_multipole) gaussian = analytic;

        EwaldVqDetail::TailSample sample;
        sample.key = key;
        sample.distance = distance;
        if (bare != nullptr) sample.bare_norm = bare->norm(2.0);
        if (gaussian != nullptr) sample.gaussian_norm = gaussian->norm(2.0);

        if (plan.coverage == TailBlockCoverage::Common && bare != nullptr && gaussian != nullptr)
        {
            RI::Tensor<Tdata> difference = *bare - *gaussian;
            sample.difference_norm = difference.norm(2.0);
            differences.emplace(key, std::move(difference));
        }
        else if (plan.coverage == TailBlockCoverage::AnalyticCancellation)
        {
            sample.bare_norm = 0.0;
            sample.gaussian_norm = 0.0;
        }
        else
        {
            local_stats.coverage_complete = false;
        }
        EwaldVqDetail::accumulate_tail_sample(local_stats, sample);
    }

    for (std::map<TailKey, TailBlockCoverage>::const_iterator coverage_it = coverage.begin();
         coverage_it != coverage.end();
         ++coverage_it)
    {
        const TailKey& key = coverage_it->first;
        const TailBlockCoverage block_coverage = coverage_it->second;
        if (!(EwaldVqDetail::canonical_hermitian_key(key) == key)) continue;
        const TailKey partner = EwaldVqDetail::hermitian_partner(key);
        const auto partner_coverage = coverage.find(partner);
        if (partner_coverage == coverage.end()
            || block_coverage == TailBlockCoverage::Incomplete
            || partner_coverage->second == TailBlockCoverage::Incomplete)
        {
            local_stats.coverage_complete = false;
            continue;
        }
        const auto difference = differences.find(key);
        const auto partner_difference = differences.find(partner);
        if (difference != differences.end() && partner_difference != differences.end())
        {
            const RI::Tensor<Tdata> partner_adjoint
                = EwaldVqDetail::tensor_adjoint(partner_difference->second);
            local_stats.hermitian_residual
                = std::max(local_stats.hermitian_residual,
                           (difference->second - partner_adjoint).norm(2.0));
        }
        else if (difference != differences.end())
        {
            local_stats.hermitian_residual
                = std::max(local_stats.hermitian_residual, difference->second.norm(2.0));
        }
        else if (partner_difference != differences.end())
        {
            local_stats.hermitian_residual
                = std::max(local_stats.hermitian_residual, partner_difference->second.norm(2.0));
        }
    }
    return EwaldVqDetail::reduce_tail_stats(this->mpi_comm, local_stats);
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_short_range_Vs(const UnitCell& ucell,
                                         const std::vector<TA>& list_A0,
                                         const std::vector<TAC>& list_A1,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_short_range_Vs");
    ModuleBase::timer::tick("Ewald_Vq", "cal_short_range_Vs");

    auto Vs_minus_gauss = this->cal_Vs_minus_gauss(ucell,
                                                   list_A0,
                                                   list_A1,
                                                   Vs_in,
                                                   this->tail_mode == EwaldVqDetail::TailMode::Enlarge
                                                       || this->tail_mode == EwaldVqDetail::TailMode::Strict);
    auto Vq_minus_gauss = this->cal_Vq_minus_gauss(ucell, list_A0, list_A1, Vs_minus_gauss);
    auto Vs_short = this->set_Vs_dVs<RI::Tensor<Tdata>>(ucell,
                                                        this->list_A0_pair_R_period,
                                                        this->list_A1_pair_R_period,
                                                        Vq_minus_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_short_range_Vs");
    return Vs_short;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_short_range_Vs_serial_full(
    const UnitCell& ucell,
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in_full,
    const std::array<Tcell, Ndim>& period_Vs_NAO) -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    const MPI_Comm mpi_comm_saved = this->mpi_comm;
    const auto list_A0_saved = this->list_A0;
    const auto list_A1_saved = this->list_A1;
    const auto list_A0_k_saved = this->list_A0_k;
    const auto list_A1_k_saved = this->list_A1_k;
    const auto list_A0_pair_R_saved = this->list_A0_pair_R;
    const auto list_A1_pair_R_saved = this->list_A1_pair_R;
    const auto list_A0_pair_R_period_saved = this->list_A0_pair_R_period;
    const auto list_A1_pair_R_period_saved = this->list_A1_pair_R_period;
    const auto list_A0_pair_k_saved = this->list_A0_pair_k;
    const auto list_A1_pair_k_saved = this->list_A1_pair_k;

    const std::array<Tcell, Ndim>& period_Vs = period_Vs_NAO;
    const std::array<int, 1> Nks = {this->nks0};

    this->mpi_comm = MPI_COMM_SELF;
    this->list_A0 = this->atoms_vec;
    this->list_A1 = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs);
    this->list_A0_k = this->atoms_vec;
    this->list_A1_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);
    this->list_A0_pair_R = this->atoms_vec;
    this->list_A1_pair_R = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs_NAO);
    this->list_A0_pair_R_period = this->atoms_vec;
    this->list_A1_pair_R_period = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, this->nmp);
    this->list_A0_pair_k = this->atoms_vec;
    this->list_A1_pair_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);

    auto Vs_short_full = this->cal_short_range_Vs(ucell, this->list_A0, this->list_A1, Vs_in_full);

    this->mpi_comm = mpi_comm_saved;
    this->list_A0 = list_A0_saved;
    this->list_A1 = list_A1_saved;
    this->list_A0_k = list_A0_k_saved;
    this->list_A1_k = list_A1_k_saved;
    this->list_A0_pair_R = list_A0_pair_R_saved;
    this->list_A1_pair_R = list_A1_pair_R_saved;
    this->list_A0_pair_R_period = list_A0_pair_R_period_saved;
    this->list_A1_pair_R_period = list_A1_pair_R_period_saved;
    this->list_A0_pair_k = list_A0_pair_k_saved;
    this->list_A1_pair_k = list_A1_pair_k_saved;

    return Vs_short_full;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_long_range_Vs_gauss(const UnitCell& ucell, const double& chi)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_long_range_Vs_gauss");
    ModuleBase::timer::tick("Ewald_Vq", "cal_long_range_Vs_gauss");

    const int shift_for_mpi = std::floor(this->nks0 / 2.0);
    T_func_DPget_Vq_dVq<RI::Tensor<std::complex<double>>> func_DPget_Vq;
    if (this->ewald_dimension == 2)
    {
        func_DPget_Vq = std::bind(&Gaussian_Abfs::get_Vq_2d,
                                  &this->gaussian_abfs,
                                  std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3,
                                  chi,
                                  std::placeholders::_4,
                                  this->gaunt);
    }
    else
    {
        func_DPget_Vq = std::bind(&Gaussian_Abfs::get_Vq,
                                  &this->gaussian_abfs,
                                  std::placeholders::_1,
                                  std::placeholders::_2,
                                  std::placeholders::_3,
                                  chi,
                                  std::placeholders::_4,
                                  this->gaunt);
    }
    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq_gauss_out
        = this->set_Vq_dVq_gauss(ucell, this->list_A0_k, this->list_A1_k, shift_for_mpi, func_DPget_Vq);

    std::set<TA> atoms10;
    std::set<TA> atoms11;
    for (const auto& I: this->list_A0_pair_k)
    {
        atoms10.insert(I);
    }
    for (const auto& JR: this->list_A1_pair_k)
    {
        atoms11.insert(JR.first);
    }
    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_gauss_out, atoms10, atoms11);

    using Tin_convert = typename LRI_CV_Tools::TinType<RI::Tensor<std::complex<double>>>::type;
    std::map<TA, std::map<TAK, RI::Tensor<std::complex<double>>>> Vq_long_full;
#pragma omp parallel
    for (size_t i0 = 0; i0 < this->list_A0_pair_k.size(); ++i0)
    {
#pragma omp for schedule(dynamic) nowait
        for (size_t i1 = 0; i1 < this->list_A1_pair_k.size(); ++i1)
        {
            const TA iat0 = this->list_A0_pair_k[i0];
            const int it0 = ucell.iat2it[iat0];
            const TA iat1 = this->list_A1_pair_k[i1].first;
            const int it1 = ucell.iat2it[iat1];
            const int ik = this->list_A1_pair_k[i1].second[0] + shift_for_mpi;
            const TAK re_index = std::make_pair(iat1, std::array<int, 1>{ik});

            auto it_outer = Vq_gauss.find(this->list_A0_pair_k[i0]);
            if (it_outer == Vq_gauss.end())
            {
                continue;
            }
            auto it_inner = it_outer->second.find(this->list_A1_pair_k[i1]);
            if (it_inner == it_outer->second.end())
            {
                continue;
            }

            RI::Tensor<std::complex<double>> data;
            LRI_CV_Tools::init_elem(data, this->index_abfs[it0].count_size, this->index_abfs[it1].count_size);
            for (int l0 = 0; l0 != this->g_abfs_ccp.at(it0).size(); ++l0)
            {
                for (int l1 = 0; l1 != this->g_abfs.at(it1).size(); ++l1)
                {
                    for (size_t n0 = 0; n0 != this->g_abfs_ccp.at(it0).at(l0).size(); ++n0)
                    {
                        const double pA = this->multipole.at(it0).at(l0).at(n0);
                        for (size_t n1 = 0; n1 != this->g_abfs.at(it1).at(l1).size(); ++n1)
                        {
                            const double pB = this->multipole.at(it1).at(l1).at(n1);
                            Tin_convert frac = RI::Global_Func::convert<Tin_convert>(pA * pB);
                            for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                            {
                                const size_t index0 = this->index_abfs.at(it0).at(l0).at(n0).at(m0);
                                const size_t lm0 = l0 * l0 + m0;
                                for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                {
                                    const size_t index1 = this->index_abfs.at(it1).at(l1).at(n1).at(m1);
                                    const size_t lm1 = l1 * l1 + m1;
                                    LRI_CV_Tools::add_elem(data,
                                                           index0,
                                                           index1,
                                                           it_inner->second,
                                                           lm0,
                                                           lm1,
                                                           frac);
                                }
                            }
                        }
                    }
                }
            }

#pragma omp critical(Ewald_Vq_cal_long_range_Vs_gauss)
            Vq_long_full[this->list_A0_pair_k[i0]][re_index] = std::move(data);
        }
    }

    auto Vs_long = this->set_Vs_dVs<RI::Tensor<Tdata>>(ucell,
                                                       this->list_A0_pair_R_period,
                                                       this->list_A1_pair_R_period,
                                                       Vq_long_full);

    ModuleBase::timer::tick("Ewald_Vq", "cal_long_range_Vs_gauss");
    return Vs_long;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_long_range_Vs_gauss_serial_full(const UnitCell& ucell,
                                                          const double& chi,
                                                          const std::array<Tcell, Ndim>& period_Vs_NAO)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    const MPI_Comm mpi_comm_saved = this->mpi_comm;
    const auto list_A0_saved = this->list_A0;
    const auto list_A1_saved = this->list_A1;
    const auto list_A0_k_saved = this->list_A0_k;
    const auto list_A1_k_saved = this->list_A1_k;
    const auto list_A0_pair_R_saved = this->list_A0_pair_R;
    const auto list_A1_pair_R_saved = this->list_A1_pair_R;
    const auto list_A0_pair_R_period_saved = this->list_A0_pair_R_period;
    const auto list_A1_pair_R_period_saved = this->list_A1_pair_R_period;
    const auto list_A0_pair_k_saved = this->list_A0_pair_k;
    const auto list_A1_pair_k_saved = this->list_A1_pair_k;

    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times, ucell, this->g_lcaos_rcut);
    const std::array<int, 1> Nks = {this->nks0};

    this->mpi_comm = MPI_COMM_SELF;
    this->list_A0 = this->atoms_vec;
    this->list_A1 = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs);
    this->list_A0_k = this->atoms_vec;
    this->list_A1_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);
    this->list_A0_pair_R = this->atoms_vec;
    this->list_A1_pair_R = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, period_Vs_NAO);
    this->list_A0_pair_R_period = this->atoms_vec;
    this->list_A1_pair_R_period = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, this->nmp);
    this->list_A0_pair_k = this->atoms_vec;
    this->list_A1_pair_k = RI::Divide_Atoms::traversal_atom_period(this->atoms_vec, Nks);

    auto Vs_long_full = this->cal_long_range_Vs_gauss(ucell, chi);

    this->mpi_comm = mpi_comm_saved;
    this->list_A0 = list_A0_saved;
    this->list_A1 = list_A1_saved;
    this->list_A0_k = list_A0_k_saved;
    this->list_A1_k = list_A1_k_saved;
    this->list_A0_pair_R = list_A0_pair_R_saved;
    this->list_A1_pair_R = list_A1_pair_R_saved;
    this->list_A0_pair_R_period = list_A0_pair_R_period_saved;
    this->list_A1_pair_R_period = list_A1_pair_R_period_saved;
    this->list_A0_pair_k = list_A0_pair_k_saved;
    this->list_A1_pair_k = list_A1_pair_k_saved;

    return Vs_long_full;
}

template <typename Tdata>
template <typename Tout, typename Tin>
auto Ewald_Vq<Tdata>::set_Vs_dVs(const UnitCell& ucell,
                                 const std::vector<TA>& list_A0_pair_R,
                                 const std::vector<TAC>& list_A1_pair_R,
                                 std::map<TA, std::map<TAK, Tin>>& Vq) -> std::map<TA, std::map<TAC, Tout>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vs_dVs");
    ModuleBase::timer::start("Ewald_Vq", "set_Vs_dVs");

    using namespace RI::Array_Operator;
    using Tin_convert = typename LRI_CV_Tools::TinType<Tout>::type;

    const double cfrac = 1.0 / this->nks0;
    std::map<TA, std::map<TAC, Tout>> datas;

    // Accumulate one (iat0, iat1, R) block at a time so OpenMP threads do
    // not each retain a near-complete output map on dense k meshes.
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (size_t i0 = 0; i0 < list_A0_pair_R.size(); ++i0)
    {
        for (size_t i1 = 0; i1 < list_A1_pair_R.size(); ++i1)
        {
            const TA iat0 = list_A0_pair_R[i0];
            const auto it_outer = Vq.find(iat0);
            if (it_outer == Vq.end())
            {
                continue;
            }

            const TA iat1 = list_A1_pair_R[i1].first;
            const TC& cell1 = list_A1_pair_R[i1].second;

            bool has_value = false;
            Tout data;

            for (size_t ik = 0; ik != this->nks0; ++ik)
            {
                const std::complex<double> frac
                    = std::exp(-ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT
                               * (this->kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1) * ucell.latvec)))
                      * cfrac;

                const TAK index = std::make_pair(iat1, std::array<int, 1>{static_cast<int>(ik)});

                // check the Fourier transformed V(q)
                // whether ccp_rmesh_times * Rcut >= rIJ
                // skip some IJ pairs
                const auto it_inner = it_outer->second.find(index);
                if (it_inner == it_outer->second.end())
                {
                    continue;
                }

                if (LRI_CV_Tools::exist(it_inner->second))
                {
                    Tout Vq_tmp = LRI_CV_Tools::convert<Tin_convert>(LRI_CV_Tools::mul2(frac, it_inner->second));
                    if (!has_value)
                    {
                        data = std::move(Vq_tmp);
                        has_value = true;
                    }
                    else
                    {
                        data = data + Vq_tmp;
                    }
                }
            }

            if (has_value)
            {
#pragma omp critical(Ewald_Vq_set_Vs_dVs)
                { datas[iat0][list_A1_pair_R[i1]] = std::move(data); }
            }
        }
    }
    // auto end = std::chrono::system_clock::now();
    // auto duration
    //     = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "set_Vs_dVs Time: "
    //           << double(duration.count())
    //                  * std::chrono::microseconds::period::num
    //                  / std::chrono::microseconds::period::den
    //           << " s" << std::endl;

    ModuleBase::timer::end("Ewald_Vq", "set_Vs_dVs");
    return datas;
}

template <typename Tdata>
std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Ewald_Vq<Tdata>::init_gauss(
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& orb_in)
{
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> gauss;
    gauss.resize(orb_in.size());
    for (size_t T = 0; T != orb_in.size(); ++T)
    {
        gauss[T].resize(orb_in[T].size());
        for (size_t L = 0; L != orb_in[T].size(); ++L)
        {
            gauss[T][L].resize(orb_in[T][L].size());
            for (size_t N = 0; N != orb_in[T][L].size(); ++N)
            {
                gauss[T][L][N] = this->gaussian_abfs.Gauss(orb_in[T][L][N], this->ewald_lambda);
            }
        }
    }

    return gauss;
}

#endif

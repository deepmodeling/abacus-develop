//=======================
// AUTHOR : jiyy
// DATE :   2024-03-08
//=======================

#ifndef EWALD_VQ_HPP
#define EWALD_VQ_HPP

#include <RI/comm/mix/Communicate_Tensors_Map_Judge.h>
#include <RI/distribute/Distribute_Equally.h>
#include <RI/global/Global_Func-1.h>

// #include <chrono>
#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "conv_coulomb_pot_k.h"
#include "exx_abfs-abfs_index.h"
#include "exx_abfs-construct_orbs.h"
#include "gaussian_abfs.h"
#include "module_base/element_basis_index.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "singular_value.h"

#include <cmath>

template <typename Tdata>
void Ewald_Vq<Tdata>::init(const UnitCell& ucell,
                           const LCAO_Orbitals& orb,
                           const MPI_Comm& mpi_comm_in,
                           const K_Vectors* kv_in,
                           std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& lcaos_in,
                           std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in,
                           const std::map<std::string, double>& parameter,
                           ORB_gaunt_table& MGT_in)
{
    ModuleBase::TITLE("Ewald_Vq", "init");
    ModuleBase::timer::tick("Ewald_Vq", "init");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = kv_in;
    this->nks0 = this->p_kv->get_nkstot_full() / this->nspin0;
    this->kvec_c.resize(this->nks0);

    this->parameter = parameter;
    this->g_lcaos = this->init_gauss(lcaos_in);
    this->g_abfs = this->init_gauss(abfs_in);
    this->g_abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->g_abfs,
                                                        this->info.ccp_type,
                                                        this->parameter,
                                                        this->info.ccp_rmesh_times);
    this->multipole = Exx_Abfs::Construct_Orbs::get_multipole(abfs_in);
    this->lcaos_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(lcaos_in);
    this->g_lcaos_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(this->g_lcaos);
    this->g_abfs_ccp_rcut = Exx_Abfs::Construct_Orbs::get_Rcut(this->g_abfs_ccp);

    const ModuleBase::Element_Basis_Index::Range range_abfs = Exx_Abfs::Abfs_Index::construct_range(abfs_in);
    this->index_abfs = ModuleBase::Element_Basis_Index::construct_index(range_abfs);

    this->cv
        .set_orbitals(ucell, orb, this->g_lcaos, this->g_abfs, this->g_abfs_ccp, this->info.kmesh_times, MGT_in, false, false);
    this->gaunt.create(MGT_in.Gaunt_Coefficients.getBound1(),
                       MGT_in.Gaunt_Coefficients.getBound2(),
                       MGT_in.Gaunt_Coefficients.getBound3());
    this->gaunt = MGT_in.Gaunt_Coefficients;

    this->atoms_vec.resize(ucell.nat);
    std::iota(this->atoms_vec.begin(), this->atoms_vec.end(), 0);
    this->nmp = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

    ModuleBase::timer::tick("Ewald_Vq", "init");
}

template <typename Tdata>
void Ewald_Vq<Tdata>::init_ions(const UnitCell& ucell, const std::array<Tcell, Ndim>& period_Vs_NAO)
{
    ModuleBase::TITLE("Ewald_Vq", "init_ions");
    ModuleBase::timer::tick("Ewald_Vq", "init_ions");

    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, this->g_lcaos_rcut);

    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, this->atoms_vec, period_Vs, 2, false);

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
    this->gaussian_abfs.init(2 * GlobalC::exx_info.info_ri.abfs_Lmax + 1, neg_kvec, ucell.G, this->ewald_lambda);

    ModuleBase::timer::tick("Ewald_Vq", "init_ions");
}

template <typename Tdata>
double Ewald_Vq<Tdata>::get_singular_chi(const UnitCell& ucell, const Singular_Value::Fq_type& fq_type, const double& qdiv)
{
    ModuleBase::TITLE("Ewald_Vq", "get_singular_chi");
    ModuleBase::timer::tick("Ewald_Vq", "get_singular_chi");

    double chi = 0.0;
    switch (fq_type)
    {
    case Singular_Value::Fq_type::Type_0:
        chi = Singular_Value::cal_type_0(ucell, this->kvec_c, qdiv, 100, 30, 1e-6, 3);
        break;
    case Singular_Value::Fq_type::Type_1:
        chi = Singular_Value::cal_type_1(ucell, this->nmp, qdiv, 1, 5, 1e-4);
        break;
    default:
        throw std::domain_error(std::string(__FILE__) + " line " + std::to_string(__LINE__));
        break;
    }

    ModuleBase::timer::tick("Ewald_Vq", "get_singular_chi");
    return chi;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_gauss(const UnitCell& ucell, const std::vector<TA>& list_A0, const std::vector<TAC>& list_A1)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vs_gauss");
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_gauss");

    std::map<std::string, bool> flags = {{"writable_Vws", true}};
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_gauss = this->cv.cal_Vs(ucell, list_A0, list_A1, flags);
    this->cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_gauss");
    return Vs_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVs_gauss(const UnitCell& ucell, const std::vector<TA>& list_A0, const std::vector<TAC>& list_A1)
    -> std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVs_gauss");
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs_gauss");

    std::map<std::string, bool> flags = {{"writable_dVws", true}};

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_gauss
        = this->cv.cal_dVs(list_A0, list_A1, flags);
    this->cv.dVws = LRI_CV_Tools::get_dCVws(dVs_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs_gauss");
    return dVs_gauss;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_Vs_minus_gauss(const UnitCell& ucell,
                                         const std::vector<TA>& list_A0,
                                         const std::vector<TAC>& list_A1,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_in)
    -> std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_Vs_minus_gauss");
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_minus_gauss");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_gauss = this->cal_Vs_gauss(ucell, list_A0, list_A1);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_minus_gauss
        = this->set_Vs_dVs_minus_gauss(ucell, list_A0, list_A1, Vs_in, Vs_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vs_minus_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs_minus_gauss");

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_gauss = this->cal_dVs_gauss(list_A0, list_A1);
    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_minus_gauss
        = this->set_Vs_dVs_minus_gauss(ucell, list_A0, list_A1, dVs_in, dVs_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs_minus_gauss");
    return dVs_minus_gauss;
}

template <typename Tdata>
double Ewald_Vq<Tdata>::cal_V_Rcut(const int it0, const int it1)
{
    return this->g_abfs_ccp_rcut[it0] + this->g_lcaos_rcut[it1];
}

template <typename Tdata>
double Ewald_Vq<Tdata>::get_Rcut_max(const int it0, const int it1)
{
    double lcaos_rmax = this->lcaos_rcut[it0] * this->info.ccp_rmesh_times + this->lcaos_rcut[it1];
    double g_lcaos_rmax = this->g_lcaos_rcut[it0] * this->info.ccp_rmesh_times + this->g_lcaos_rcut[it1];
    return std::min(lcaos_rmax, g_lcaos_rmax);
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

    using Tin_convert = typename LRI_CV_Tools::TinType<Tresult>::type;
    std::map<TA, std::map<TAC, Tresult>> pVs_dVs_gauss;
#pragma omp parallel
    for (size_t i0 = 0; i0 < list_A0.size(); ++i0)
    {
#pragma omp for schedule(dynamic) nowait
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
                for (int l0 = 0; l0 != this->g_abfs_ccp[it0].size(); ++l0)
                {
                    for (int l1 = 0; l1 != this->g_abfs[it1].size(); ++l1)
                    {
                        for (size_t n0 = 0; n0 != this->g_abfs_ccp[it0][l0].size(); ++n0)
                        {
                            const double pA = this->multipole[it0][l0][n0];
                            for (size_t n1 = 0; n1 != this->g_abfs[it1][l1].size(); ++n1)
                            {
                                const double pB = this->multipole[it1][l1][n1];
                                Tin_convert pp = RI::Global_Func::convert<Tin_convert>(pA * pB);
                                for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                                {
                                    for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                    {
                                        const size_t index0 = this->index_abfs[it0][l0][n0][m0];
                                        const size_t index1 = this->index_abfs[it1][l1][n1][m1];

                                        LRI_CV_Tools::add_elem(data,
                                                               index0,
                                                               index1,
                                                               Vs_dVs_gauss_in[list_A0[i0]][list_A1[i1]],
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

    std::map<TA, std::map<TAC, Tresult>> Vs_dVs_minus_gauss = LRI_CV_Tools::minus(Vs_dVs_in, pVs_dVs_gauss);
    ModuleBase::timer::tick("Ewald_Vq", "set_Vs_dVs_minus_gauss");
    return Vs_dVs_minus_gauss;
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
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq_gauss");

    const T_func_DPget_Vq_dVq<RI::Tensor<std::complex<double>>> func_DPget_Vq = std::bind(&Gaussian_Abfs::get_Vq,
                                                                                          &this->gaussian_abfs,
                                                                                          ucell,
                                                                                          std::placeholders::_1,
                                                                                          std::placeholders::_2,
                                                                                          std::placeholders::_3,
                                                                                          chi,
                                                                                          std::placeholders::_4,
                                                                                          this->gaunt);
    auto Vq_gauss = this->set_Vq_dVq_gauss(ucell, list_A0_k, list_A1_k, shift_for_mpi, func_DPget_Vq);

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq_gauss");

    const T_func_DPget_Vq_dVq<std::array<RI::Tensor<std::complex<double>>, Ndim>> func_DPget_dVq
        = std::bind(&Gaussian_Abfs::get_dVq,
                    &this->gaussian_abfs,
                    ucell,
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
                dVq_gauss[iat0][index0] = -res[list_A0_k[i0]][list_A1_k[i1]];
                const TAK index1 = std::make_pair(iat0, TK{ik});
                dVq_gauss[iat1][index1] = res[list_A0_k[i0]][list_A1_k[i1]];
            }
            else
                dVq_gauss[list_A0_k[i0]][list_A1_k[i1]] = res[list_A0_k[i0]][list_A1_k[i1]];
        }
    }

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq_gauss");

    std::map<TA, std::map<TAK, Tresult>> Vq_dVq_gauss_out;

#pragma omp parallel
    for (size_t i0 = 0; i0 < list_A0_k.size(); ++i0)
    {
#pragma omp for schedule(dynamic) nowait
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
            auto data = func_DPget_Vq_dVq(this->g_abfs_ccp[it0].size() - 1, this->g_abfs[it1].size() - 1, ik, tau);

#pragma omp critical(Ewald_Vq_set_Vq_dVq_gauss)
            Vq_dVq_gauss_out[list_A0_k[i0]][list_A1_k[i1]] = data;
        }
    }

    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq_minus_gauss");

    auto Vq_minus_gauss
        = this->set_Vq_dVq_minus_gauss<RI::Tensor<std::complex<double>>>(ucell, list_A0, list_A1, Vs_minus_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq_minus_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq_minus_gauss");

    auto dVq_minus_gauss
        = this->set_Vq_dVq_minus_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>>(ucell,
                                                                                           list_A0,
                                                                                           list_A1,
                                                                                           dVs_minus_gauss);

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq_minus_gauss");
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
    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq_minus_gauss");

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
                    const double Rcut = std::min(this->get_Rcut_max(it0, it1), this->get_Rcut_max(it1, it0));
                    const ModuleBase::Vector3<double> R_delta
                        = -tau0 + tau1 + (RI_Util::array3_to_Vector3(cell1) * ucell.latvec);

                    if (R_delta.norm() * ucell.lat0 < Rcut)
                    {
                        const std::complex<double> phase
                            = std::exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT
                                       * (this->kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1) * ucell.latvec)));

                        Tout Vs_dVs_tmp = LRI_CV_Tools::mul2(
                            phase,
                            LRI_CV_Tools::convert<Tin_convert>(std::move(Vs_dVs_minus_gauss[iat0][list_A1[i1]])));

                        const TAK index = std::make_pair(iat1, TK{static_cast<int>(ik)});
                        if (!LRI_CV_Tools::exist(local_datas[iat0][index]))
                            local_datas[iat0][index] = Vs_dVs_tmp;
                        else
                            local_datas[iat0][index] = local_datas[iat0][index] + Vs_dVs_tmp;
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
                        datas[key0][key1] = datas[key0][key1] + value;
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
    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq_minus_gauss");
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
                                                                                             Vs_in); //{ia0, {ia1, R}}
    const T_func_DPcal_Vq_dVq_minus_gauss<RI::Tensor<std::complex<double>>, RI::Tensor<Tdata>> func_cal_Vq_minus_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_Vq_minus_gauss,
                    this,
                    ucell,
                    this->list_A0_pair_R,
                    this->list_A1_pair_R,
                    std::placeholders::_1);
    const T_func_DPcal_Vq_dVq_gauss<RI::Tensor<std::complex<double>>> func_cal_Vq_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_Vq_gauss,
                    this,
                    ucell,
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

    ModuleBase::timer::tick("Ewald_Vq", "cal_Vq");
    return Vq;
}

template <typename Tdata>
auto Ewald_Vq<Tdata>::cal_dVq(const UnitCell& ucell,
                              const double& chi,
                              std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_in)
    -> std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVq");
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq");

    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs_minus_gauss
        = this->cal_dVs_minus_gauss(ucell,
                                    this->list_A0,
                                    this->list_A1,
                                    dVs_in); //{ia0, {ia1, R}}
    const T_func_DPcal_Vq_dVq_minus_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>,
                                          std::array<RI::Tensor<Tdata>, Ndim>>
        func_cal_dVq_minus_gauss = std::bind(&Ewald_Vq<Tdata>::cal_dVq_minus_gauss,
                                             this,
                                             ucell,
                                             this->list_A0_pair_R,
                                             this->list_A1_pair_R,
                                             std::placeholders::_1);
    const T_func_DPcal_Vq_dVq_gauss<std::array<RI::Tensor<std::complex<double>>, Ndim>> func_cal_dVq_gauss
        = std::bind(&Ewald_Vq<Tdata>::cal_dVq_gauss,
                    this,
                    ucell,
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

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVq");
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
    for (const auto& outerPair: Vs_dVs_minus_gauss_in)
    {
        atoms00.insert(outerPair.first);

        for (const auto& innerPair: outerPair.second)
            atoms01.insert(innerPair.first.first);
    }
    std::map<TA, std::map<TAC, Tin>> Vs_dVs_minus_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_dVs_minus_gauss_in, atoms00, atoms01);
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_minus_gauss = func_cal_Vq_dVq_minus_gauss(Vs_dVs_minus_gauss); //{ia0, ia1}

    // MPI: {ia0, {ia1, k}} to {ia0, ia1}
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_gauss_out = func_cal_Vq_dVq_gauss(shift_for_mpi); //{ia0, {ia1, k}}
    std::set<TA> atoms10;
    std::set<TA> atoms11;
    for (const auto& outerPair: Vq_dVq_gauss_out)
    {
        atoms10.insert(outerPair.first);

        for (const auto& innerPair: outerPair.second)
            atoms11.insert(innerPair.first.first);
    }
    std::map<TA, std::map<TAK, Tout>> Vq_dVq_gauss
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vq_dVq_gauss_out, atoms10, atoms11); //{ia0, ia1}

#pragma omp parallel
    for (size_t i0 = 0; i0 < list_A0_pair_k.size(); ++i0)
    {
#pragma omp for schedule(dynamic) nowait
        for (size_t i1 = 0; i1 < list_A1_pair_k.size(); ++i1)
        {
            const TA iat0 = list_A0_pair_k[i0];
            const int it0 = ucell.iat2it[iat0];
            const TA iat1 = list_A1_pair_k[i1].first;
            const int it1 = ucell.iat2it[iat1];
            const int ik = list_A1_pair_k[i1].second[0] + shift_for_mpi;
            const TAK re_index = std::make_pair(iat1, std::array<int, 1>{ik});

            const size_t size0 = this->index_abfs[it0].count_size;
            const size_t size1 = this->index_abfs[it1].count_size;
            Tout data;
            LRI_CV_Tools::init_elem(data, size0, size1);
            for (int l0 = 0; l0 != this->g_abfs_ccp[it0].size(); ++l0)
            {
                for (int l1 = 0; l1 != this->g_abfs[it1].size(); ++l1)
                {
                    for (size_t n0 = 0; n0 != this->g_abfs_ccp[it0][l0].size(); ++n0)
                    {
                        const double pA = this->multipole[it0][l0][n0];
                        for (size_t n1 = 0; n1 != this->g_abfs[it1][l1].size(); ++n1)
                        {
                            const double pB = this->multipole[it1][l1][n1];
                            Tin_convert pp = RI::Global_Func::convert<Tin_convert>(pA * pB);
                            for (size_t m0 = 0; m0 != 2 * l0 + 1; ++m0)
                            {
                                const size_t index0 = this->index_abfs[it0][l0][n0][m0];
                                const size_t lm0 = l0 * l0 + m0;
                                for (size_t m1 = 0; m1 != 2 * l1 + 1; ++m1)
                                {
                                    const size_t index1 = this->index_abfs[it1][l1][n1][m1];
                                    const size_t lm1 = l1 * l1 + m1;

                                    LRI_CV_Tools::add_elem(data,
                                                           index0,
                                                           index1,
                                                           Vq_dVq_gauss[list_A0_pair_k[i0]][list_A1_pair_k[i1]],
                                                           lm0,
                                                           lm1,
                                                           pp);
                                }
                            }
                        }
                    }
                }
            }

#pragma omp critical(Ewald_Vq_set_Vq_dVq)
            if (LRI_CV_Tools::exist(Vq_dVq_minus_gauss[list_A0_pair_k[i0]][re_index]))
                Vq_dVq[list_A0_pair_k[i0]][re_index] = Vq_dVq_minus_gauss[list_A0_pair_k[i0]][re_index] + data;
        }
    }

    ModuleBase::timer::tick("Ewald_Vq", "set_Vq_dVq");
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
auto Ewald_Vq<Tdata>::cal_dVs(
    const UnitCell& ucell,
    const double& chi,
    std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>& dVs_in) //{ia0, {ia1, R}}
    -> std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
{
    ModuleBase::TITLE("Ewald_Vq", "cal_dVs");
    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs");

    std::map<TA, std::map<TAK, std::array<RI::Tensor<std::complex<double>>, Ndim>>> dVq
        = this->cal_dVq(ucell, chi, dVs_in);
    auto dVs = this->set_Vs_dVs<std::array<RI::Tensor<Tdata>, Ndim>>(ucell,
                                                                     this->list_A0_pair_R_period,
                                                                     this->list_A1_pair_R_period,
                                                                     dVq); //{ia0, ia1}

    ModuleBase::timer::tick("Ewald_Vq", "cal_dVs");
    return dVs;
}

template <typename Tdata>
template <typename Tout, typename Tin>
auto Ewald_Vq<Tdata>::set_Vs_dVs(const UnitCell& ucell,
                                 const std::vector<TA>& list_A0_pair_R,
                                 const std::vector<TAC>& list_A1_pair_R,
                                 std::map<TA, std::map<TAK, Tin>>& Vq) -> std::map<TA, std::map<TAC, Tout>>
{
    ModuleBase::TITLE("Ewald_Vq", "set_Vs_dVs");
    ModuleBase::timer::tick("Ewald_Vq", "set_Vs_dVs");

    using namespace RI::Array_Operator;
    using Tin_convert = typename LRI_CV_Tools::TinType<Tout>::type;

    const double cfrac = 1.0 / this->nks0;
    std::map<TA, std::map<TAC, Tout>> datas;

    // auto start = std::chrono::system_clock::now();
#pragma omp parallel
    {
        std::map<TA, std::map<TAC, Tout>> local_datas;

#pragma omp for schedule(dynamic) nowait
        for (size_t ik = 0; ik != this->nks0; ++ik)
        {
            for (size_t i0 = 0; i0 < list_A0_pair_R.size(); ++i0)
            {
                for (size_t i1 = 0; i1 < list_A1_pair_R.size(); ++i1)
                {
                    const TA iat0 = list_A0_pair_R[i0];
                    const TA iat1 = list_A1_pair_R[i1].first;
                    const TC& cell1 = list_A1_pair_R[i1].second;
                    const std::complex<double> frac
                        = std::exp(-ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT
                                   * (this->kvec_c[ik] * (RI_Util::array3_to_Vector3(cell1) * ucell.latvec)))
                          * cfrac;

                    const TAK index = std::make_pair(iat1, std::array<int, 1>{static_cast<int>(ik)});
                    if (LRI_CV_Tools::exist(Vq[iat0][index]))
                    {
                        Tout Vq_tmp = LRI_CV_Tools::convert<Tin_convert>(LRI_CV_Tools::mul2(frac, Vq[iat0][index]));

                        if (!LRI_CV_Tools::exist(local_datas[iat0][list_A1_pair_R[i1]]))
                            local_datas[iat0][list_A1_pair_R[i1]] = Vq_tmp;
                        else
                            local_datas[iat0][list_A1_pair_R[i1]] = local_datas[iat0][list_A1_pair_R[i1]] + Vq_tmp;
                    }
                }
            }
        }
#pragma omp critical(Ewald_Vq_set_Vs_dVs)
        {
            for (auto& outer_pair: local_datas)
            {
                TA key0 = outer_pair.first;
                for (auto& inner_pair: outer_pair.second)
                {
                    TAC key1 = inner_pair.first;
                    Tout& value = inner_pair.second;
                    if (!LRI_CV_Tools::exist(datas[key0][key1]))
                        datas[key0][key1] = value;
                    else
                        datas[key0][key1] = datas[key0][key1] + value;
                }
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

    ModuleBase::timer::tick("Ewald_Vq", "set_Vs_dVs");
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

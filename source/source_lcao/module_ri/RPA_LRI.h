//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_H
#define RPA_LRI_H

#include "source_esolver/esolver_ks_lcao.h"
#include "LRI_CV.h"
// #include "module_xc/exx_info.h"
// #include "source_basis/module_ao/ORB_atomic_lm.h"
#include "source_base/matrix.h"
// #include "source_lcao/module_ri/Exx_LRI.h"
// #include <RI/physics/Exx.h>
#include <RI/ri/RI_Tools.h>
#include <array>
#include <map>
#include <mpi.h>
#include <vector>

class Parallel_Orbitals;
class K_Vectors;

template <typename T, typename Tdata> class RPA_LRI
{
  private:
    using TA = int;
    using Tcell = int;
    static constexpr std::size_t Ndim = 3;
    using TC = std::array<Tcell, Ndim>;
    using Tq = std::array<double, Ndim>;
    using TAC = std::pair<TA, TC>;
    using TAq = std::pair<TA, Tq>;
    using TatomR = std::array<double, Ndim>; // tmp

  public:
    RPA_LRI(const Exx_Info::Exx_Info_RI &info_in) : info(info_in)
    {
    }
    ~RPA_LRI(){};
    void postSCF(const UnitCell& ucell,
        const MPI_Comm& mpi_comm_in,
        const elecstate::DensityMatrix<T, Tdata>& dm,
        const elecstate::ElecState* pelec,
        const K_Vectors& kv,
        const LCAO_Orbitals& orb,
        const Parallel_Orbitals& parav,
        const psi::Psi<T>& psi);
    void init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in, const std::vector<double>& orb_cutoff);
    void cal_postSCF_exx(const elecstate::DensityMatrix<T, Tdata>& dm,
        const MPI_Comm& mpi_comm_in,
        const UnitCell& ucell,
        const K_Vectors& kv,
        const LCAO_Orbitals& orb);
    void output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb);
    void cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv);
    void cal_abfs_overlap(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv);
    void inverse_olp(const UnitCell& ucell,
                     std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                     const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
    void out_abfs_overlap(const UnitCell& ucell,
                          std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                          std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                          std::string filename,
                          const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                          const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
    void out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi);
    void out_struc(const UnitCell& ucell);
    void out_bands(const elecstate::ElecState *pelec);

    void output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa);
    void out_Cs(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename);
    void out_coulomb_k(const UnitCell& ucell,
                       std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                       std::string filename,
                       Exx_LRI<double>* exx_lri);
    // void print_matrix(char *desc, const ModuleBase::matrix &mat);
    // void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);
    // void init(const MPI_Comm &mpi_comm_in);
    // void cal_rpa_ions();

    Tdata Erpa;

  private:
    const Exx_Info::Exx_Info_RI &info;
    const K_Vectors *p_kv=nullptr;
    MPI_Comm mpi_comm;
    std::vector<double> orb_cutoff_;
    double ccp_rmesh_times_ewald = 0.0;

    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
    // shrinked abfs
    std::shared_ptr<ORB_gaunt_table> MGT;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_shrink;


    // Exx_LRI<double> exx_postSCF_double(info);
    // LRI_CV<Tdata> cv;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_period;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs_period;
    // RI::RPA<TA,Tcell,Ndim,Tdata> rpa_lri;

    // Tdata post_process_Erpa( const Tdata &Erpa_in ) const;

    Exx_LRI<double>* exx_cut_coulomb = nullptr;
    Exx_LRI<double>* exx_full_coulomb = nullptr;
};
Exx_LRI<double> exx_lri_rpa(GlobalC::exx_info.info_ri);
#include "RPA_LRI.hpp"

#endif

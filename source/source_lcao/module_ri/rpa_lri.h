//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_H
#define RPA_LRI_H

#include "LRI_CV.h"
#include "source_esolver/esolver_ks_lcao.h"
#include "source_lcao/module_ri/module_exx_symmetry/symmetry_rotation.h"
#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"
// #include "module_xc/exx_info.h"
// #include "source_basis/module_ao/orb_atomic_lm.h"
#include "source_base/matrix.h"
// #include "source_lcao/module_ri/exx_lri.h"
// #include <RI/physics/Exx.h>
#include <RI/ri/RI_Tools.h>
#include <array>
#include <map>
#include <memory>
#include <mpi.h>
#include <vector>

class Parallel_Orbitals;
class K_Vectors;
template <typename Tdata> class Exx_LRI;

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
    RPA_LRI(const Exx_Info_RI &info_in) : info(info_in)
    {
    }
    ~RPA_LRI();
    void postSCF(const UnitCell& ucell,
        const MPI_Comm& mpi_comm_in,
        const module_dm::DensityMatrix<T, Tdata>& dm,
        const elecstate::ElecState* pelec,
        const K_Vectors& kv,
        const LCAO_Orbitals& orb,
        const Parallel_Orbitals& parav,
        const psi::Psi<T>& psi);
    void init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in, const std::vector<double>& orb_cutoff);
    void cal_postSCF_exx(const module_dm::DensityMatrix<T, Tdata>& dm,
        const MPI_Comm& mpi_comm_in,
        const UnitCell& ucell,
        const K_Vectors& kv,
        const LCAO_Orbitals& orb,
        const Parallel_Orbitals& parav);
    void output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb);
    void cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv);
    void cal_abfs_overlap(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv);
    void output_symmetry_sidecars(const UnitCell& ucell,
                                  const K_Vectors& kv,
                                  const elecstate::DensityMatrix<T, Tdata>& dm);
    void inverse_olp(const UnitCell& ucell,
                     std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                     const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
    void out_abfs_overlap(const UnitCell& ucell,
                          std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                          std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                          std::string filename,
                          const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                          const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
    void out_abfs_overlap_v1(const UnitCell& ucell,
                             std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                             std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                             std::string filename,
                             const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                             const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
    void out_abfs_overlap_raw_v1(const UnitCell& ucell,
                                 const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                 const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
    void out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi);
    void out_struc(const UnitCell& ucell);
    void out_bz_sampling();
    void out_bands(const elecstate::ElecState *pelec);

    ModuleRI::SternheimerOrbitalSet take_sternheimer_abfs();
    static void trim_process_heap();

    void output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa);
    void out_Cs(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename);
    void out_Cs_v1(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename);
    void out_coulomb_k(const UnitCell& ucell,
                       std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                       std::string filename,
                       Exx_LRI<double>* exx_lri);
    void out_coulomb_k_v1(const UnitCell& ucell,
                          std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                          std::string filename,
                          Exx_LRI<double>* exx_lri);
    void out_librpa_basis_v1(const UnitCell& ucell,
                             Exx_LRI<double>* exx_lri,
                             const std::string& aux_filename = "basis_aux_out",
                             const std::string& legacy_filename = "basis_out");
    // void print_matrix(char *desc, const ModuleBase::matrix &mat);
    // void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);
    // void init(const MPI_Comm &mpi_comm_in);
    // void cal_rpa_ions();

    Tdata Erpa;

  private:
    Conv_Coulomb_Pot_K::Coulomb_Method select_coulomb_basis_method_(Exx_LRI<double>* exx_lri) const;
    std::vector<int> collect_atom_naux_(const UnitCell& ucell, Exx_LRI<double>* exx_lri) const;

    const Exx_Info::Exx_Info_RI &info;
    const K_Vectors *p_kv=nullptr;
    MPI_Comm mpi_comm;
    std::vector<double> orb_cutoff_;
    double ccp_rmesh_times_ewald = 0.0;
    double ccp_rmesh_times_cut = 0.0;

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

    ModuleSymmetry::Symmetry_rotation symmetry_rotation_;
    bool use_spacegroup_symmetry_ = false;
    std::unique_ptr<Exx_LRI<double>> exx_cut_coulomb;
    std::unique_ptr<Exx_LRI<double>> exx_full_coulomb;
};
#include "rpa_lri.hpp"

#endif

//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_H
#define RPA_LRI_H

#include "module_esolver/esolver_ks_lcao.h"
#include "LRI_CV.h"
// #include "module_xc/exx_info.h"
// #include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_base/matrix.h"
// #include "module_ri/Exx_LRI.h"
// #include <RI/physics/Exx.h>
#include <RI/ri/RI_Tools.h>
#include <array>
#include <map>
#include <mpi.h>
#include <vector>

class Local_Orbital_Charge;
class Parallel_Orbitals;
class K_Vectors;

template <typename T, typename Tdata> class RPA_LRI
{
  private:
    using TA = int;
    using Tcell = int;
    static constexpr std::size_t Ndim = 3;
    using TC = std::array<Tcell, Ndim>;
    using TAC = std::pair<TA, TC>;
    using TatomR = std::array<double, Ndim>; // tmp

  public:
    RPA_LRI(const Exx_Info::Exx_Info_RI &info_in) : info(info_in)
    {
    }
    ~RPA_LRI(){};
    void init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in);
    void cal_rpa_cv();
    void cal_postSCF_exx(const elecstate::DensityMatrix<T, Tdata>& dm,
        const MPI_Comm& mpi_comm_in,
        const K_Vectors& kv);
    void out_for_RPA(const Parallel_Orbitals& parav,
        const psi::Psi<T>& psi,
        const elecstate::ElecState* pelec);
    void out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi);
    void out_struc();
    void out_bands(const elecstate::ElecState *pelec);

    void out_Cs();
    void out_coulomb_k();
    // void print_matrix(char *desc, const ModuleBase::matrix &mat);
    // void print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat);
    // void init(const MPI_Comm &mpi_comm_in);
    // void cal_rpa_ions();
    // void cal_rpa_elec(const Local_Orbital_Charge &loc, const Parallel_Orbitals &pv);

    Tdata Erpa;

  private:
    const Exx_Info::Exx_Info_RI &info;
    const K_Vectors *p_kv=nullptr;
    MPI_Comm mpi_comm;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;

    // Exx_LRI<double> exx_postSCF_double(info);
    // LRI_CV<Tdata> cv;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_period;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs_period;
    // RI::RPA<TA,Tcell,Ndim,Tdata> rpa_lri;

    // Tdata post_process_Erpa( const Tdata &Erpa_in ) const;
};
Exx_LRI<double> exx_lri_rpa(GlobalC::exx_info.info_ri);
#include "RPA_LRI.hpp"

#endif
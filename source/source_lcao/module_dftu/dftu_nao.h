#ifndef DFTU_LCAO_H
#define DFTU_LCAO_H

#include "source_pw/module_pwdft/dftu_base.h"

#include <complex>
#include <string>
#include <vector>


class UnitCell;
class Parallel_Orbitals;

#ifdef __LCAO
class LCAO_Orbitals;
namespace hamilt
{
template <typename T>
class HContainer;
} // namespace hamilt
namespace elecstate
{
template <typename TK, typename TR>
class DensityMatrix;
} // namespace elecstate
#endif


class Plus_U : public Plus_U_Base
{

  public:
    Plus_U();
    ~Plus_U();

    // allocate relevant data strcutures
    void init(UnitCell& cell,
                const Parallel_Orbitals* pv,
                const int npol,
                const int nspin,
                const std::vector<int>& orbital_corr,
                const bool yukawa_potential,
                const double yukawa_lambda,
                const std::string& global_readin_dir,
                const std::string& global_out_dir,
                const std::string& init_chg,
                const int nlocal,
                const std::string& ks_solver,
                const std::string& device,
                const int kpar,
                const std::vector<double>& hubbard_u,
                const double uramping,
                const int occ_mat_ctrl,
                const int mixing_dftu
#ifdef __LCAO
                , const LCAO_Orbitals* orb = nullptr
#endif
                );

  private:

#ifdef __LCAO
    const LCAO_Orbitals* ptr_orb_ = nullptr;
    std::vector<double> orb_cutoff_;

    //=============================================================
    // In dftu_hamilt.cpp (DFTU_LCAO free functions)
    // For calculating contribution to Hamiltonian matrices
    //=============================================================
  public:
    /**
     * @brief get the density matrix of target spin
     * nspin = 1 and 4 : ispin should be 0
     * nspin = 2 : ispin should be 0/1
    */
    const hamilt::HContainer<double>* get_dmr(int ispin) const;
    /**
     * @brief set the density matrix for DFT+U calculation
     * if the density matrix is not set or set to nullptr, the DFT+U calculation will not be performed
    */
    void set_dmr(const elecstate::DensityMatrix<double, double>* dm_in_dftu_d);
    void set_dmr(const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd);

    /// read-only accessors for state needed by DFTU_LCAO free functions
    const std::vector<double>& get_orb_cutoff() const { return orb_cutoff_; }
    const LCAO_Orbitals* get_ptr_orb() const { return ptr_orb_; }

  private:
    const elecstate::DensityMatrix<double, double>* dm_in_dftu_d = nullptr;
    const elecstate::DensityMatrix<std::complex<double>, double>* dm_in_dftu_cd = nullptr;
#endif
};

#endif

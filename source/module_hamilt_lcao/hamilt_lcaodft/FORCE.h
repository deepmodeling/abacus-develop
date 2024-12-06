#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_FORCE_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HAMILT_LCAO_HAMILT_LCAODFT_FORCE_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/force_stress_arrays.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_psi/psi.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_elecstate/elecstate.h"

#ifndef TGINT_H
#define TGINT_H
template <typename T>
struct TGint;
template <>
struct TGint<double>
{
    using type = Gint_Gamma;
};
template <>
struct TGint<std::complex<double>>
{
    using type = Gint_k;
};
#endif

template <typename T>
class Force_Stress_LCAO;

template <typename T>
class Force_LCAO
{
  public:
    friend class Force_Stress_LCAO<T>;

    Force_LCAO(){};
    ~Force_LCAO(){};

  private:
    const Parallel_Orbitals* ParaV;

    elecstate::Potential* pot;

    // orthonormal force + contribution from T and VNL
    void ftable(const bool isforce,
                const bool isstress,
                ForceStressArrays& fsr, // mohan add 2024-06-16
                const UnitCell& ucell,
                const psi::Psi<T>* psi,
                const elecstate::ElecState* pelec,
                ModuleBase::matrix& foverlap,
                ModuleBase::matrix& ftvnl_dphi,
                ModuleBase::matrix& fvnl_dbeta,
                ModuleBase::matrix& fvl_dphi,
                ModuleBase::matrix& soverlap,
                ModuleBase::matrix& stvnl_dphi,
                ModuleBase::matrix& svnl_dbeta,
                ModuleBase::matrix& svl_dphi,
#ifdef __DEEPKS
                ModuleBase::matrix& svnl_dalpha,
#endif
                typename TGint<T>::type& gint,
                const TwoCenterBundle& two_center_bundle,
                const LCAO_Orbitals& orb,
                const Parallel_Orbitals& pv,
                const K_Vectors* kv = nullptr,
                Record_adj* ra = nullptr);

    // get the ds, dt, dvnl.
    void allocate(const UnitCell& ucell,
                  const Parallel_Orbitals& pv,
                  ForceStressArrays& fsr, // mohan add 2024-06-15
                  const TwoCenterBundle& two_center_bundle,
                  const LCAO_Orbitals& orb,
                  const int& nks = 0,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_d = {});

    void finish_ftable(ForceStressArrays& fsr);

    void average_force(double* fm);

    //void test(Parallel_Orbitals& pv, double* mm, const std::string& name);

    //-------------------------------------------------------------
    // forces reated to overlap matrix
    // forces related to energy density matrix
    //-------------------------------------------------------------

    void cal_fedm(const bool isforce,
                  const bool isstress,
                  ForceStressArrays& fsr,
                  const UnitCell& ucell,
                  const elecstate::DensityMatrix<T, double>& dm,
                  const psi::Psi<T>* psi,
                  const Parallel_Orbitals& pv,
                  const elecstate::ElecState* pelec,
                  ModuleBase::matrix& foverlap,
                  ModuleBase::matrix& soverlap,
                  const K_Vectors* kv = nullptr,
                  Record_adj* ra = nullptr);

    //-------------------------------------------------------------
    // forces related to kinetic and non-local pseudopotentials
    //--------------------------------------------------------------
    void cal_ftvnl_dphi(const elecstate::DensityMatrix<T, double>* dm,
                        const Parallel_Orbitals& pv,
                        const UnitCell& ucell,
                        ForceStressArrays& fsr,
                        const bool isforce,
                        const bool isstress,
                        ModuleBase::matrix& ftvnl_dphi,
                        ModuleBase::matrix& stvnl_dphi,
                        Record_adj* ra = nullptr);

    //-------------------------------------------
    // forces related to local pseudopotentials
    //-------------------------------------------
    void cal_fvl_dphi(const bool isforce,
                      const bool isstress,
                      const elecstate::Potential* pot_in,
                      typename TGint<T>::type& gint,
                      ModuleBase::matrix& fvl_dphi,
                      ModuleBase::matrix& svl_dphi);

    elecstate::DensityMatrix<T, double> cal_edm(const elecstate::ElecState* pelec,
        const psi::Psi<T>& psi,
        const elecstate::DensityMatrix<T, double>& dm,
        const K_Vectors& kv,
        const Parallel_Orbitals& pv,
        const int& nspin, 
        const int& nbands,
        const UnitCell& ucell,
        Record_adj& ra) const;
};

#endif

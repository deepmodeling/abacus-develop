#ifndef FORCE_LCAO_GAMMA_H
#define FORCE_LCAO_GAMMA_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_psi/psi.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"


#ifndef TGINT_H
#define TGINT_H
template <typename T> struct TGint;
template <>
struct TGint<double> {
    using type = Gint_Gamma;
};
template <>
struct TGint<std::complex<double>> {
    using type = Gint_k;
};
#endif

template<typename T> class Force_Stress_LCAO;

template<typename T>
class Force_LCAO
{
public:
    friend class Force_Stress_LCAO<T>;

    Force_LCAO() {};
    ~Force_LCAO() {};

private:

    const Parallel_Orbitals* ParaV;

    elecstate::Potential* pot;

    // orthonormal force + contribution from T and VNL
    void ftable(
        const bool isforce,
        const bool isstress,
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
        const ORB_gen_tables* uot,
        const Parallel_Orbitals& pv,
        LCAO_Matrix& lm,
        const K_Vectors* kv = nullptr,
        Record_adj* ra = nullptr);


    // get the ds, dt, dvnl.
    void allocate(const Parallel_Orbitals& pv,
        LCAO_Matrix& lm,
        const ORB_gen_tables* uot,
        const int& nks = 0,
        const std::vector<ModuleBase::Vector3<double>>& kvec_d = {});


    void finish_ftable(LCAO_Matrix& lm);

    void average_force(double* fm);

    void test(Parallel_Orbitals& pv, double* mm, const std::string& name);

    //-------------------------------------------------------------
    // forces reated to overlap matrix
    // forces related to energy density matrix
    //-------------------------------------------------------------

    void cal_fedm(
        const bool isforce,
        const bool isstress,
        const UnitCell& ucell,
        const elecstate::DensityMatrix<T, double>* dm,
        const psi::Psi<T>* psi,
        const Parallel_Orbitals& pv,
        const elecstate::ElecState* pelec,
        LCAO_Matrix& lm,
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
        LCAO_Matrix& lm,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix& ftvnl_dphi,
        ModuleBase::matrix& stvnl_dphi,
        Record_adj* ra = nullptr);

    void cal_fvnl_dbeta(const elecstate::DensityMatrix<T, double>* dm,
        const Parallel_Orbitals& pv,
        const UnitCell& ucell,
        const LCAO_Orbitals& orb,
        const ORB_gen_tables& uot,
        Grid_Driver& gd,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix& fvnl_dbeta,
        ModuleBase::matrix& svnl_dbeta);

    //-------------------------------------------
    // forces related to local pseudopotentials
    //-------------------------------------------
    void cal_fvl_dphi(const bool isforce,
        const bool isstress,
        const elecstate::Potential* pot_in,
        typename TGint<T>::type& gint,
        ModuleBase::matrix& fvl_dphi,
        ModuleBase::matrix& svl_dphi);
};

// this namespace used to store global function for some stress operation
namespace StressTools
{
    // set upper matrix to whole matrix
    void stress_fill(const double& lat0_, const double& omega_, ModuleBase::matrix& stress_matrix);
} // namespace StressTools
#endif

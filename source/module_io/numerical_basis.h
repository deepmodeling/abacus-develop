//==========================================================
// AUTHOR : mohan
// DATE : 2009-4-2
// Last Modify: 2009-08-28
//==========================================================
#ifndef NUMERICAL_BASIS_H
#define NUMERICAL_BASIS_H
#include <vector>

#include "bessel_basis.h"
#include "module_base/complexarray.h"
#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/vector3.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
//==========================================================
// CLASS :
// NAME :  Numerical_Basis
//==========================================================
class Numerical_Basis
{
  public:
    // this is not a good idea to construct the instance. Instead, there are many variables that CAN DEFINE the identity of this instance,
    // these parameters should be provided in the constructor. For future refactor, I list them here:
    // bessel_nao_*: including ecut, rcut, smearing, sigma, etc.
    // a file name: for the function start_from_file_k, if starts from file, can construct a different instance, instead of using the same instance.
    Numerical_Basis();
    ~Numerical_Basis();

    // void start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf, const ModulePW::PW_Basis_K* wfcpw, const UnitCell& ucell);
    void output_overlap(const psi::Psi<std::complex<double>>& psi,
                        const Structure_Factor& sf,
                        const K_Vectors& kv,
                        const ModulePW::PW_Basis_K* wfcpw,
                        const UnitCell& ucell,
                        const int& index);

  private:
    bool init_label = false;

    Bessel_Basis bessel_basis;

    std::vector<ModuleBase::IntArray> mu_index;
    static std::vector<ModuleBase::IntArray> init_mu_index(const UnitCell& ucell);

    void numerical_atomic_wfc(const int& ik,
                              const ModulePW::PW_Basis_K* wfcpw,
                              ModuleBase::ComplexMatrix& psi,
                              const Structure_Factor& sf,
                              const UnitCell& ucell);

    ModuleBase::ComplexArray cal_overlap_Q(const int& ik,
                                           const int& np,
                                           const ModulePW::PW_Basis_K* wfcpw,
                                           const psi::Psi<std::complex<double>>& psi,
                                           const double derivative_order,
                                           const Structure_Factor& sf,
                                           const UnitCell& ucell) const;

    // computed in the plane-wave basis
    ModuleBase::ComplexArray cal_overlap_Sq(const int& ik,
                                            const int& np,
                                            const double derivative_order,
                                            const Structure_Factor& sf,
                                            const ModulePW::PW_Basis_K* wfcpw,
                                            const UnitCell& ucell) const;

    static ModuleBase::matrix cal_overlap_V(const ModulePW::PW_Basis_K* wfcpw,
                                            const psi::Psi<std::complex<double>>& psi,
                                            const double derivative_order,
                                            const K_Vectors& kv,
                                            const double tpiba);

    // gk should be in the atomic unit (Bohr)
    ModuleBase::realArray cal_flq(const std::vector<ModuleBase::Vector3<double>> &gk,
                                  const int ucell_lmax) const;

    // Ylm does not depend on the magnitude so unit is not important
    static ModuleBase::matrix cal_ylm(const std::vector<ModuleBase::Vector3<double>> &gk,
                                      const int ucell_lmax);

    // gk and the returned gpow are both in the atomic unit (Bohr)
    static std::vector<double> cal_gpow(const std::vector<ModuleBase::Vector3<double>> &gk,
                                        const double derivative_order);

    static void output_info(std::ofstream& ofs, const Bessel_Basis& bessel_basis, const K_Vectors& kv, const UnitCell& ucell);

    static void output_k(std::ofstream& ofs, const K_Vectors& kv);

    static void output_overlap_Q(std::ofstream& ofs,
                                 const std::vector<ModuleBase::ComplexArray>& overlap_Q,
                                 const K_Vectors& kv);

    static void output_overlap_Sq(const std::string& name,
                                  std::ofstream& ofs,
                                  const std::vector<ModuleBase::ComplexArray>& overlap_Sq,
                                  const K_Vectors& kv);

    static void output_overlap_V(std::ofstream &ofs, const ModuleBase::matrix &overlap_V);
};

#endif

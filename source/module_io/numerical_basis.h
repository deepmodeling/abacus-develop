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
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
//==========================================================
// CLASS :
// NAME :  Numerical_Basis
//==========================================================
class Numerical_Basis
{
  public:
    Numerical_Basis();
    ~Numerical_Basis();

    void start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor* p_sf);
    void output_overlap(const psi::Psi<std::complex<double>>& psi, const Structure_Factor* p_sf);

  private:
    bool init_label = false;

    Bessel_Basis bessel_basis;

    std::vector<ModuleBase::IntArray> mu_index;
    static std::vector<ModuleBase::IntArray> init_mu_index(void);

    void numerical_atomic_wfc(const int& ik,
                              ModulePW::PW_Basis_K* wfc_basis,
                              ModuleBase::ComplexMatrix& psi,
                              const Structure_Factor* p_sf);

    ModuleBase::ComplexArray cal_overlap_Q(const int& ik,
                                           const int& np,
                                           ModulePW::PW_Basis_K* wfc_basis,
                                           const psi::Psi<std::complex<double>>& psi,
                                           const double derivative_order,
                                           const Structure_Factor* p_sf) const;

    ModuleBase::ComplexArray cal_overlap_Sq(const int& ik,
                                            const int& np,
                                            const double derivative_order,
                                            const Structure_Factor* p_sf) const;

    static ModuleBase::matrix cal_overlap_V(ModulePW::PW_Basis_K *wfc_basis,
                                            const psi::Psi<std::complex<double>> &psi,
                                            const double derivative_order);

    ModuleBase::realArray cal_flq(const int ik, const std::vector<ModuleBase::Vector3<double>> &gk) const;

    static ModuleBase::matrix cal_ylm(const std::vector<ModuleBase::Vector3<double>> &gk);

    static std::vector<double> cal_gpow(const std::vector<ModuleBase::Vector3<double>> &gk,
                                        const double derivative_order);

    static void output_info(std::ofstream &ofs, const Bessel_Basis &bessel_basis);

    static void output_k(std::ofstream &ofs);

    static void output_overlap_Q(std::ofstream &ofs, const std::vector<ModuleBase::ComplexArray> &overlap_Q);

    static void output_overlap_Sq(const std::string &name,
                                  std::ofstream &ofs,
                                  const std::vector<ModuleBase::ComplexArray> &overlap_Sq);

    static void output_overlap_V(std::ofstream &ofs, const ModuleBase::matrix &overlap_V);
};

#endif

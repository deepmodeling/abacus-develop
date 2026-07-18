//==========================================================
// AUTHOR : mohan
// DATE : 2009-4-2
// Last Modify: 2009-08-28
//==========================================================
#ifndef NUMERICAL_BASIS_H
#define NUMERICAL_BASIS_H
#include <complex>
#include <string>
#include <vector>

#include "bessel_basis.h"
#include "../../source_base/complexarray.h"
#include "../../source_base/complexmatrix.h"
#include "../../source_base/global_function.h"
#include "../../source_base/intarray.h"
#include "../../source_base/matrix.h"
#include "../../source_base/vector3.h"
#include "../../source_basis/module_pw/pw_basis_k.h"
#include "../../source_cell/klist.h"
#include "../../source_pw/module_pwdft/stru_fac.h"
#include "../../source_psi/psi.h"
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

    struct SIABPrimitiveParameters
    {
        double ecut_ry;
        double rcut_bohr;
        bool smooth;
        double sigma;
        double tolerance;
    };

    struct SIABPrimitiveGridBlock
    {
        int type_index;
        std::string element;
        int atom_index;
        int l;
        int m;
        int n_primitive;
        int offset;
        /// Rank-local values ordered by primitive index: values[ie][local_grid_index].
        std::vector<std::vector<std::complex<double>>> values;
    };

    /// Return the ABACUS real-spherical-harmonic index for a conventional m.
    static int siab_abacus_m(const int conventional_m);

    /// Read the numerical-Bessel settings used by the existing spillage path.
    static SIABPrimitiveParameters siab_parameters_from_input(const int rcut_index);

    /// Build rank-local, physically normalized SIAB primitives on the PW FFT grid.
    std::vector<SIABPrimitiveGridBlock> siab_primitive_grid_values(
        const int ik,
        const ModulePW::PW_Basis_K* wfcpw,
        const Structure_Factor& sf,
        const UnitCell& ucell,
        const SIABPrimitiveParameters& parameters);

    // void start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf, const ModulePW::PW_Basis_K* wfcpw, const UnitCell& ucell);
    void output_overlap(const psi::Psi<std::complex<double>>& psi,
                        const Structure_Factor& sf,
                        const K_Vectors& kv,
                        const ModulePW::PW_Basis_K* wfcpw,
                        const UnitCell& ucell,
                        const int& index);

  private:
    bool init_label = false;
    SIABPrimitiveParameters initialized_siab_parameters = {0.0, 0.0, false, 0.0, 0.0};
    int initialized_siab_ntype = 0;
    int initialized_siab_lmax = 0;
    int initialized_siab_nmax = 0;
    std::vector<int> initialized_siab_basis_shape;

    Bessel_Basis bessel_basis;

    std::vector<ModuleBase::IntArray> mu_index;
    static std::vector<ModuleBase::IntArray> init_mu_index(const UnitCell& ucell);
    static std::vector<int> siab_basis_shape_signature(const UnitCell& ucell);

    void initialize_siab_basis(const UnitCell& ucell, const SIABPrimitiveParameters& parameters);

    static void fill_reciprocal_primitive(
        const int l,
        const int abacus_m,
        const int primitive_index,
        const std::complex<double>* structure_factor,
        const ModuleBase::realArray& flq,
        const ModuleBase::matrix& ylm,
        const std::vector<double>& gpow,
        const double normalization,
        std::vector<std::complex<double>>& values);

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

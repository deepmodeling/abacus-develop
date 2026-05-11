#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <complex>
#include <map>
#include <vector>

#include "source_base/constants.h"
#include "source_base/complexmatrix.h"
#include "source_base/matrix.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/vector3.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/operator.h"
#include "source_estate/elecstate.h"

#ifdef __LCAO
#include "source_estate/module_dm/density_matrix.h" // mohan add 2025-11-02
#endif

namespace spinconstrain
{

/**
 * @brief Extract magnetic moment from nspin=4 occupation matrix elements.
 *
 * Given occ[4] = {|a|^2, a* b, b* a, |b|^2} (spinor density matrix),
 * the magnetic moment components are:
 *   Mz = occ[0] - occ[3]  (sigma_z)
 *   Mx = occ[1] + occ[2]  (sigma_x)
 *   My = Im(occ[1] - occ[2])  (sigma_y)
 */
inline ModuleBase::Vector3<double> pauli_to_moment(const std::complex<double> occ[4], double weight)
{
    return ModuleBase::Vector3<double>(
        weight * (occ[1] + occ[2]).real(),
        weight * (occ[1] - occ[2]).imag(),
        weight * (occ[0] - occ[3]).real()
    );
}

struct ScAtomData;

template <typename TK>
class SpinConstrain
{
public:
    /**
     * pubic interface for spin-constrained DFT
    */
    /// initialize spin-constrained DFT
  void init_sc(double sc_thr_in,
               int nsc_in,
               int nsc_min_in,
               double alpha_trial_in,
               double sccut_in,
               double sc_drop_thr_in,
               const UnitCell& ucell,
               bool direction_only_in,
               Parallel_Orbitals* ParaV_in,
               int nspin_in,
               const K_Vectors& kv_in,
               void* p_hamilt_in,
               void* psi_in,
#ifdef __LCAO
			   elecstate::DensityMatrix<TK, double> *dm_in, // mohan add 2025-11-02
#endif
			   elecstate::ElecState* pelec_in,
               ModulePW::PW_Basis_K* pw_wfc_in = nullptr);

  /// @brief calculate the magnetization of each atom with real space projection method for LCAO base
  /// @param step : the step number of the SCF calculation
  /// @param print : print the magnetization of each atom if true
  void cal_mi_lcao(const int& step, bool print = false);

  void cal_mi_pw();

  void cal_mw_from_lambda(int i_step, 
		  const ModuleBase::Vector3<double>* delta_lambda = nullptr);

  /**
   * @brief calculate the energy of \sum_i \lambda_i * Mi
   * if this->is_mag_converged is true, then this function will calculate the energy and return the real value
   * if this->is_mag_converged is false, then this function will return 0.0
   */
  double cal_escon();

  double get_escon() const;

  void run_lambda_loop(int outer_step,
		  bool rerun = true);

  void run_lambda_bfgs_v2(int outer_step);

  void run_lambda_linear_scan(int outer_step);

  /// @brief update the charge density for LCAO base with new lambda
  /// update the charge density and psi for PW base with new lambda
  void update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve = true, bool full_update = false);

  /**
   * @brief Wavefunction and charge density update implementation for PW basis
   * @details Two-stage process:
   *          1. Subspace diagonalization: apply DeltaSpin correction and solve for each k-point
   *          2. Charge update: full-space diagonalization or direct charge update based on pw_solve
   */
  void update_psi_charge_pw(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
  
  /// CPU implementation of PW basis update
  void update_psi_charge_pw_cpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
  
#if ((defined __CUDA) || (defined __ROCM))
  /// GPU implementation of PW basis update
  void update_psi_charge_pw_gpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
#endif

  void calculate_delta_hcc(std::complex<double>* h_tmp,
		  const std::complex<double>* becp_k,
		  const ModuleBase::Vector3<double>* delta_lambda,
		  const int nbands, const int nkb, const int* nh_iat, const int ik,
		  bool full_update = false);

#ifdef __LCAO
  /// @brief convert orbital matrix to nested vector format
  std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix& orbMulP);
  /// @brief calculate magnetic moment from orbital matrix
  void calculate_MW(const std::vector<std::vector<std::vector<double>>>& AorbMulP);
  /// @brief collect magnetic moment from complex matrix
  void collect_MW(ModuleBase::matrix& MecMulP,
                  const ModuleBase::ComplexMatrix& mud,
                  int nw,
                  int isk);
#endif

  /// lambda loop helper functions
  bool check_rms_stop(int outer_step, int i_step, double rms_error, double duration, double total_duration);

  /// apply restriction
  void check_restriction(const std::vector<ModuleBase::Vector3<double>>& search, double& alpha_trial);

  /// check gradient decay
  bool check_gradient_decay(std::vector<ModuleBase::Vector3<double>> new_spin,
                            std::vector<ModuleBase::Vector3<double>> old_spin,
                            std::vector<ModuleBase::Vector3<double>> new_delta_lambda,
                            std::vector<ModuleBase::Vector3<double>> old_delta_lambda,
                            bool print = false);
  /// @brief  calculate alpha_opt
  double cal_alpha_opt(std::vector<ModuleBase::Vector3<double>> spin,
                       std::vector<ModuleBase::Vector3<double>> spin_plus,
                       const double alpha_trial);
  /// print header info
  void print_header();
  /// print termination message
  void print_termination();

  /// print mi
  void print_Mi(std::ofstream& ofs_running);

  /// print magnetic force, defined as \frac{\delta{L}}/{\delta{Mi}} = -lambda[iat])
  void print_Mag_Force(std::ofstream& ofs_running);

  /// @brief use rerun to get higher precision in lambda_loop for PW base
  bool higher_mag_prec = false;

public:
    /**
     * important outter class pointers used in spin-constrained DFT
    */
    Parallel_Orbitals *ParaV = nullptr;
    //--------------------------------------------------------------------------------
    // pointers for solve Hamiltonian to get new Magnetization from Lambda
    void* p_hamilt = nullptr;
    void* psi = nullptr;
    elecstate::ElecState* pelec = nullptr;
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
#ifdef __LCAO
    elecstate::DensityMatrix<TK, double>* dm_;
#endif
    double tpiba = 0.0; /// save ucell.tpiba
    const double meV_to_Ry = 7.349864435130999e-05;
    K_Vectors kv_;
    //--------------------------------------------------------------------------------

  public:
    /**
     * pubic methods for setting and getting spin-constrained DFT parameters
    */
    /// Public method to access the Singleton instance
    static SpinConstrain& getScInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// set element index to atom index map
    void set_atomCounts(const std::map<int, int>& atomCounts_in);
    /// get element index to atom index map
    const std::map<int, int>& get_atomCounts() const;
    /// set element index to orbital index map
    void set_orbitalCounts(const std::map<int, int>& orbitalCounts_in);
    /// get element index to orbital index map
    const std::map<int, int>& get_orbitalCounts() const;
    /// set lnchiCounts
    void set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in);
    /// get lnchiCounts
    const std::map<int, std::map<int, int>>& get_lnchiCounts() const;
    /// set sc_lambda
    void set_sc_lambda();
    /// set sc_lambda from variable
    void set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in);
    /// set target_mag
    void set_target_mag();
    /// set target_mag from variable
    void set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in);
    /// set target magnetic moment
    void set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in);
    /// set constrain
    void set_constrain();
    /// set constrain from variable
    void set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in);
    /// get sc_lambda
    const std::vector<ModuleBase::Vector3<double>>& get_sc_lambda() const;
    /// get target_mag
    const std::vector<ModuleBase::Vector3<double>>& get_target_mag() const;
    /// get constrain
    const std::vector<ModuleBase::Vector3<int>>& get_constrain() const;
    /// get nat
    int get_nat();
    /// get ntype
    int get_ntype();
    /// check atomCounts
    void check_atomCounts();
    /// get iat
    int get_iat(int itype, int atom_index);
    /// set nspin
    void set_nspin(int nspin);
    /// get nspin
    int get_nspin() const;
    /// zero atomic magnetic moment
    void zero_Mi();
    /// get decay_grad
    double get_decay_grad(int itype);
    /// set decay_grad
    void set_decay_grad();
    /// get decay_grad
    const std::vector<double>& get_decay_grad();
    /// set decay_grad from variable
    void set_decay_grad(const double* decay_grad_in, int ntype_in);
    /// set decay grad switch
    void set_sc_drop_thr(double sc_drop_thr_in);
    /// set input parameters
    void set_input_parameters(double sc_thr_in,
                              int nsc_in,
                              int nsc_min_in,
                              double alpha_trial_in,
                              double sccut_in,
                              double sc_drop_thr_in);
    /// get sc_thr
    double get_sc_thr() const;
    /// get nsc
    int get_nsc() const;
    /// get nsc_min
    int get_nsc_min() const;
    /// get alpha_trial
    double get_alpha_trial() const;
    /// get sccut
    double get_sccut() const;
    /// get sc_drop_thr
    double get_sc_drop_thr() const;
    /// @brief set orbital parallel info
    void set_ParaV(Parallel_Orbitals* ParaV_in);
    /// @brief set parameters for solver
    void set_solver_parameters(const K_Vectors& kv_in,
                               void* p_hamilt_in,
                               void* psi_in,
                               elecstate::ElecState* pelec_in);

  private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain()
    {
        delete[] sub_h_save;
        delete[] sub_s_save;
        delete[] becp_save;
        sub_h_save = nullptr;
        sub_s_save = nullptr;
        becp_save = nullptr;
    };
    SpinConstrain& operator=(SpinConstrain const&) = delete;  // Copy assign
    SpinConstrain& operator=(SpinConstrain &&) = delete;      // Move assign
    std::map<int, std::vector<ScAtomData>> ScData;
    std::map<int, double> ScDecayGrad; // in unit of uB^2/eV
    std::vector<double> decay_grad_;   // in unit of uB^2/Ry
    std::map<int, int> atomCounts;
    std::map<int, int> orbitalCounts;
    std::map<int, std::map<int, int>> lnchiCounts;
    std::vector<ModuleBase::Vector3<double>> lambda_; // in unit of Ry/uB in code, but in unit of meV/uB in input file
    std::vector<ModuleBase::Vector3<double>> target_mag_; // in unit of uB
    std::vector<ModuleBase::Vector3<double>> Mi_; // in unit of uB
    std::vector<std::string> atomLabels_;
    double escon_ = 0.0;
    int nspin_ = 0;
    int npol_ = 1;
    /**
     * input parameters for lambda-loop
     */
    int nsc_;
    int nsc_min_;
    double sc_drop_thr_ = 1e-3;
    double sc_thr_; // in unit of uB
    double current_sc_thr_;
    std::vector<ModuleBase::Vector3<int>> constrain_;
    bool debug = false;
    double alpha_trial_; // in unit of Ry/uB^2 = 0.01 eV/uB^2
    double restrict_current_; // in unit of Ry/uB = 3 eV/uB
    bool direction_only_ = false; ///< only optimize the direction of magnetization

  public:
    /// @brief save operator for spin-constrained DFT
    /// @param op_in the base pointer of operator, actual type should be DeltaSpin<OperatorLCAO<TK, TR>>*
    void set_operator(hamilt::Operator<TK>* op_in);
    /// @brief set is_Mi_converged
    void set_mag_converged(bool is_Mi_converged_in){this->is_Mi_converged = is_Mi_converged_in;}
    /// @brief get is_Mi_converged
    bool mag_converged() const {return this->is_Mi_converged;}
    void set_npol(int npol);
    int get_npol() const;
    int get_nw() const;
    int get_iwt(int itype, int iat, int orbital_index) const;
    /// get spin sign for k-point ik: +1 for spin-up, -1 for spin-down
    int get_spin_sign(int ik) const;
    /// accumulate Mi from becp for a single k-point
    void accumulate_Mi_from_becp(const std::complex<double>* becp,
                                 int nkb,
                                 int nbands,
                                 int npol,
                                 int ik,
                                 const double* wg_ik,
                                 const int* nh_iat);
  private:
    /// operator for spin-constrained DFT, used for calculating current atomic magnetic moment
    hamilt::Operator<TK>* p_operator = nullptr;
    /// @brief if atomic magnetic moment is converged
    bool is_Mi_converged = false;

    TK* sub_h_save = nullptr;
    TK* sub_s_save = nullptr;
    TK* becp_save = nullptr;
    std::vector<ModuleBase::Vector3<double>> lambda_in_sub_;
};


/**
 * @brief struct for storing parameters of non-collinear spin-constrained DFT
 */
struct ScAtomData {
    int index;
    std::vector<double> lambda;
    std::vector<double> target_mag;
    std::vector<int> constrain;
    int mag_type;
    double target_mag_val;
    double target_mag_angle1;
    double target_mag_angle2;
};

} // namespace spinconstrain

#endif // SPIN_CONSTRAIN_H

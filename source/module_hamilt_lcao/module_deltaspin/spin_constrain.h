#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <map>
#include <vector>

#include "module_base/constants.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"
#include "module_base/vector3.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hsolver/hsolver.h"

struct ScAtomData;

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
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
               bool decay_grad_switch_in,
               const UnitCell& ucell,
               std::string sc_file,
               int NPOL,
               Parallel_Orbitals* ParaV_in,
               int nspin_in,
               K_Vectors kv_in,
               std::string KS_SOLVER_in,
               LCAO_Matrix* LM_in,
               hsolver::HSolver<FPTYPE, Device>* phsol_in,
               hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
               psi::Psi<FPTYPE>* psi_in,
               elecstate::ElecState* pelec_in);

  /// calculate h_lambda operator for spin-constrained DFT
  void cal_h_lambda(std::complex<double>* h_lambda, const std::vector<std::complex<double>>& Sloc2, bool column_major);

  void cal_MW(const int& step, LCAO_Matrix* LM, const UnitCell& ucell, bool print = false);

  ModuleBase::matrix cal_MW_k(LCAO_Matrix* LM, const std::vector<std::vector<std::complex<double>>>& dm);

  void cal_mw_from_lambda(int i_step);

  double cal_escon();

  double get_escon();

  std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix& orbMulP);

  void run_lambda_loop(int outer_step);

public:
    /**
     * important outter class pointers used in spin-constrained DFT
    */
    Parallel_Orbitals *ParaV = nullptr;
    hsolver::HSolver<FPTYPE, Device>* phsol = nullptr;
    hamilt::Hamilt<FPTYPE, Device>* p_hamilt = nullptr;
    psi::Psi<FPTYPE>* psi = nullptr;
    elecstate::ElecState* pelec = nullptr;
    LCAO_Matrix* LM = nullptr;
    std::string KS_SOLVER;
    const double meV_to_Ry = 7.349864435130999e-05;
    K_Vectors kv_;

  public:
    /**
     * pubic methods for setting and getting spin-constrained DFT parameters
    */
    /// Public method to access the Singleton instance
    static SpinConstrain& getScInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// parse json input file for non-collinear spin-constrained DFT
    void Set_ScData_From_Json(const std::string& filename);
    /// get sc_data
    const std::map<int, std::vector<ScAtomData>>& get_ScData() const;
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
    /// get nw
    int get_nw();
    /// get iwt
    int get_iwt(int itype, int iat, int orbital_index);
    /// set npol
    void set_npol(int npol);
    /// get npol
    int get_npol();
    /// set nspin
    void set_nspin(int nspin);
    /// get nspin
    int get_nspin();
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
    void set_decay_grad_switch(bool decay_grad_switch_in);
    /// set input parameters
    void set_input_parameters(double sc_thr_in,
                              int nsc_in,
                              int nsc_min_in,
                              double alpha_trial_in,
                              double sccut_in,
                              bool decay_grad_switch_in);
    /// get sc_thr
    double get_sc_thr();
    /// get nsc
    int get_nsc();
    /// get nsc_min
    int get_nsc_min();
    /// get alpha_trial
    double get_alpha_trial();
    /// get sccut
    double get_sccut();
    /// get decay_grad_switch
    bool get_decay_grad_switch();
    /// @brief set orbital parallel info
    void set_ParaV(Parallel_Orbitals* ParaV_in);
    /// @brief set parameters for solver
    void set_solver_parameters(int nspin_in,
                               K_Vectors kv_in,
                               hsolver::HSolver<FPTYPE, Device>* phsol_in,
                               hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                               psi::Psi<FPTYPE>* psi_in,
                               elecstate::ElecState* pelec_in,
                               std::string KS_SOLVER_in,
                               LCAO_Matrix* LM_in);
    /// bcast sc data read from json file
    void bcast_ScData(std::string sc_file, int nat, int ntype);

  private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain(){};                              // Destructor
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
    double escon_ = 0.0;
    int nspin_ = 0;
    int npol_ = 1;
    /**
     * input parameters for lambda-loop
     */
    int nsc_;
    int nsc_min_;
    bool decay_grad_switch_ = false;
    double sc_thr_; // in unit of uB
    std::vector<ModuleBase::Vector3<int>> constrain_;
    bool debug = false;
    double alpha_trial_; // in unit of Ry/uB^2 = 0.01 eV/uB^2
    double restrict_current_; // in unit of Ry/uB = 3 eV/uB
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

#endif // SPIN_CONSTRAIN_H

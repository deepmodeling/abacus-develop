#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <vector>
#include <map>
#include "module_base/vector3.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/unitcell.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hsolver/hsolver.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"

struct ScAtomData;

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class SpinConstrain
{
public:
    /**
     * pubic interface for spin-constrained DFT
    */
    /// initialize spin-constrained DFT
    void init_sc(const UnitCell& ucell,
                int NPOL,
                std::string sc_file,
                Parallel_Orbitals* ParaV_in,
                int nspin_in,
                double sc_thr_in,
                int nsc_in,
                int nsc_min_in,
                K_Vectors kv_in,
                std::string KS_SOLVER_in,
                LCAO_Matrix* LM_in,
                hsolver::HSolver<FPTYPE, Device>* phsol_in,
                hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                psi::Psi<FPTYPE>* psi_in,
                elecstate::ElecState* pelec_in);

    /// calculate h_lambda operator for spin-constrained DFT
    void cal_h_lambda(std::complex<double>* h_lambda);

    void cal_weight_func(const std::vector<std::complex<double>>& Sloc2);

    void cal_MW(
        const int& step,
        LCAO_Matrix& LM,
        const std::vector<ModuleBase::ComplexMatrix> &dm,
        const UnitCell& ucell,
        bool print=false);

    ModuleBase::matrix cal_MW_k(
        LCAO_Matrix& LM,
        const std::vector<ModuleBase::ComplexMatrix> &dm
    );

    void cal_mw_from_lambda(int i_step);

    double cal_escon();

    std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix &orbMulP);

    inline double output_cut(const double& result)
    {
        if(std::abs(result) < 1e-6)
        {
            return 0.0;
        }
        return result;
    }

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
    /// clear sc_data
    void clear_ScData();
    /// set element index to atom index map
    void set_atomCounts(const std::map<int, int>& atomCounts_in);
    /// get element index to atom index map
    const std::map<int, int>& get_atomCounts() const;
    /// set element index to orbital index map
    void set_orbitalCounts(const std::map<int, int>& orbitalCounts_in);
    /// get element index to orbital index map
    const std::map<int, int>& get_orbitalCounts() const;
    /// set sc_lambda
    void set_sc_lambda();
    /// set sc_lambda from variable
    void set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in);
    /// set init_mag
    void set_init_mag();
    /// set init_mag from variable
    void set_init_mag(const ModuleBase::Vector3<double>* init_mag_in, int nat_in);
    /// set sc_mag
    void set_sc_mag();
    /// set sc_mag from variable
    void set_sc_mag(const ModuleBase::Vector3<double>* sc_mag_in, int nat_in);
    /// set constrain
    void set_constrain();
    /// set constrain from variable
    void set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in);
    /// get sc_lambda
    const std::vector<ModuleBase::Vector3<double>>& get_sc_lambda() const;
    /// get init_mag
    const std::vector<ModuleBase::Vector3<double>>& get_init_mag() const;
    /// get sc_mag
    const std::vector<ModuleBase::Vector3<double>>& get_sc_mag() const;
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
    /// clear atomCounts
    void clear_atomCounts();
    /// clear orbitalCounts
    void clear_orbitalCounts();
    /// set npol
    void set_npol(int npol);
    /// get npol
    int get_npol();
    /// set nspin
    void set_nspin(int nspin);
    /// get nspin
    int get_nspin();

private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain(){};                              // Destructor
    SpinConstrain& operator=(SpinConstrain const&) = delete;  // Copy assign
    SpinConstrain& operator=(SpinConstrain &&) = delete;      // Move assign
    std::map<int, std::vector<ScAtomData>> ScData;
    std::map<int, int> atomCounts;
    std::map<int, int> orbitalCounts;
    int nspin_ = 0;
    int npol_ = 1;
    std::vector<std::complex<double>> Wi_; // same as overlap matrix
    std::vector<ModuleBase::Vector3<double>> lambda_; // in unit of Ry/uB in code, but in unit of meV/uB in input file
    std::vector<ModuleBase::Vector3<double>> sc_mag_; // in unit of uB
    std::vector<ModuleBase::Vector3<double>> Mi_; // in unit of uB
    double escon_ = 0.0;
    /**
     * parameters for lambda-loop
    */
    int nsc_;
    int nsc_min_;
    double sc_thr_; // in unit of uB
    std::vector<ModuleBase::Vector3<int>> constrain_;
    std::vector<ModuleBase::Vector3<double>> init_mag_;
    bool debug = false;
    double alpha_trial_ = 0.00073498; // in unit of Ry/uB^2 = 0.01 eV/uB^2
    double restrict_current_ = 0.22049; // in unit of Ry/uB = 3 eV/uB
    K_Vectors kv_;
};


/**
 * @brief struct for storing parameters of non-collinear spin-constrained DFT
 */
struct ScAtomData {
    int index;
    std::vector<double> lambda;
    std::vector<double> init_mag;
    std::vector<double> sc_mag;
    std::vector<int> constrain;
    double sc_spin_val;
    double sc_spin_angle1;
    double sc_spin_angle2;
};

#endif // SPIN_CONSTRAIN_H

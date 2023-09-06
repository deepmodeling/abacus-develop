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
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"

struct ScAtomData;

class SpinConstrain
{
public:
    /// Public method to access the Singleton instance
    static SpinConstrain& getInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// init
    void init_sc(const UnitCell& ucell, int NPOL, std::string sc_file, Parallel_Orbitals* ParaV_in);
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
    /// set sc_mag
    void set_sc_mag();
    /// set sc_mag from variable
    void set_sc_mag(const ModuleBase::Vector3<double>* sc_mag_in, int nat_in);
    /// get sc_lambda
    const std::vector<ModuleBase::Vector3<double>>& get_sc_lambda() const;
    /// get sc_mag
    const std::vector<ModuleBase::Vector3<double>>& get_sc_mag() const;
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
    /// calculate h_lambda operator for spin-constrained DFT
    void cal_h_lambda(const std::vector<std::complex<double>>& Sloc2, std::complex<double>* h_lambda);

public:
    Parallel_Orbitals *ParaV;
    void cal_MW(const int& step,
            LCAO_Hamilt& uhm,
            Local_Orbital_Charge& loc,
            const K_Vectors& kv,
            const UnitCell& ucell);

private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain(){};                              // Destructor
    SpinConstrain& operator=(SpinConstrain const&) = delete;  // Copy assign
    SpinConstrain& operator=(SpinConstrain &&) = delete;      // Move assign
    std::map<int, std::vector<ScAtomData>> ScData;
    std::map<int, int> atomCounts;
    std::map<int, int> orbitalCounts;
    int npol_ = 1;
    std::vector<std::complex<double>> Wi_;
    std::vector<ModuleBase::Vector3<double>> lambda_;
    std::vector<ModuleBase::Vector3<double>> sc_mag_;
};


/**
 * @brief struct for storing parameters of non-collinear spin-constrained DFT
 */
struct ScAtomData {
    int index;
    std::vector<double> lambda;
    std::vector<double> sc_mag;
    double sc_spin_val;
    double sc_spin_angle1;
    double sc_spin_angle2;
};

#endif // SPIN_CONSTRAIN_H
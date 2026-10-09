#ifndef INFONONLOCAL_H
#define INFONONLOCAL_H

#include <fstream>
#include <string>
#include <vector>

#include "../source_cell/atom_spec.h"
#include "../source_basis/module_ao/orb_nonlocal.h"
#include "../source_basis/module_ao/orb_read.h"

class InfoNonlocal
{
public:
    InfoNonlocal();
    ~InfoNonlocal();

    /// NON-LOCAL part for LCAO
private:
    std::vector<Numerical_Nonlocal> Beta; ///< nonlocal projectors (one per atom type)
    std::vector<int> nproj;               ///< number of projectors per atom type, mohan add 2010-12-19
    int nprojmax;                         ///< max number of projectors among all types, mohan add 2010-03-07
    double rcutmax_Beta;                  ///< max cutoff radius among all projectors, caoyu add 2021-05-24

public:

    const std::vector<Numerical_Nonlocal>& get_Beta() const { return Beta; }
    std::vector<Numerical_Nonlocal>& get_Beta() { return Beta; }
    const Numerical_Nonlocal& get_Beta(const int& it) const { return Beta[it]; }
    const Numerical_Nonlocal* get_Beta_data() const { return Beta.data(); }
    Numerical_Nonlocal* get_Beta_data() { return Beta.data(); }
    // Replaces the whole vector without copying or moving any element,
    // since Numerical_Nonlocal is non-copyable (see orb_nonlocal.h).
    void resize_Beta(const int& ntype) { Beta = std::vector<Numerical_Nonlocal>(ntype); }

    const std::vector<int>& get_nproj() const { return nproj; }
    std::vector<int>& get_nproj() { return nproj; }
    int get_nproj(const int& it) const { return nproj[it]; }
    int& get_nproj_ref(const int& it) { return nproj[it]; }
    void assign_nproj(const int& ntype, const int& value) { nproj.assign(ntype, value); }

    const int& get_nprojmax() const { return nprojmax; }
    void set_nprojmax(const int& value) { nprojmax = value; }

    const double& get_rcutmax_Beta(void) const { return rcutmax_Beta; }
    void set_rcutmax_Beta(const double& value) { rcutmax_Beta = value; }

    /// in order to get rid of the .NONLOCAL file.
    void Set_NonLocal(
        const int& it,
        Atom* atom,
        int& n_projectors,
        const int& kmesh,
        const double& dk,
        const double& dr_uniform,
        std::ofstream& log,
        const bool& out_element_info,
        const bool& lspinorb,
        const int& nspin,
        const int& my_rank);

    /// workflow to setup nonlocal part for LCAO
    void setupNonlocal(
        const int& ntype,
        Atom* atoms,
        std::ofstream& log,
        LCAO_Orbitals& orb,
        const std::string& basis_type,
        const bool& out_element_info,
        const bool& lspinorb,
        const int& nspin,
        const int& my_rank);

private:
    /// build SOC coefficient matrix for nonlocal projectors
    void build_soc_coefficients(
        const Atom* atom,
        const int& n_projectors,
        ModuleBase::ComplexMatrix& coefficient_D_nc_in);

    /// build radial projector beta_r on truncated mesh
    void build_beta_r(
        const Atom* atom,
        const int& p1,
        std::vector<double>& beta_r,
        int& cut_mesh);
};

#endif

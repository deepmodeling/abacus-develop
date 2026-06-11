// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PHON_H
#define DFPT_PHON_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"

namespace ModuleDFPT {

class DFPT_Phon {
public:
    DFPT_Phon();
    ~DFPT_Phon();
    
    void init(UnitCell& ucell);
    
    void assemble(int q_idx, DFPT_PW_Data& data);
    
    void diagonalize(int q_idx, DFPT_PW_Data& data);
    
    void add_loto(DFPT_PW_Data& data);
    
    bool check_sum_rule(int q_idx, DFPT_PW_Data& data) const;

private:
    UnitCell* ucell_ = nullptr;
    
    double ewald_alpha_ = 0.0;
    double ewald_rcut_ = 0.0;
    
    void ion_ion(const ModuleBase::Vector3<double>& q, ModuleBase::matrix& dyn);
    
    void electron(int q_idx, DFPT_PW_Data& data, ModuleBase::matrix& dyn);
    
    void ewald_sum(const ModuleBase::Vector3<double>& q, ModuleBase::matrix& dyn);
};

} // namespace ModuleDFPT

#endif // DFPT_PHON_H
// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef QLIST_H
#define QLIST_H

#include "source_base/vector3.h"
#include "module_symmetry/symmetry.h"
#include "unitcell.h"
#include <vector>

namespace ModuleCell {

class QList {
public:
    QList();
    ~QList();
    
    void generate_mesh(UnitCell& ucell, ModuleSymmetry::Symmetry& symm,
                       const std::vector<int>& mp_grid, bool use_irreps);
    
    void read_from_file(const std::string& filename, UnitCell& ucell);
    
    int get_nq() const { return nq_; }
    
    ModuleBase::Vector3<double> get_q(int idx) const { return qvec_[idx]; }
    
    int get_nirr(int idx) const { return nirr_[idx]; }
    
    std::vector<int> get_irrep_modes(int q_idx, int irrep_idx) const;

private:
    int nq_ = 0;
    std::vector<ModuleBase::Vector3<double>> qvec_;
    std::vector<int> nirr_;
    std::vector<std::vector<std::vector<int>>> irrep_modes_;
    
    void reduce(UnitCell& ucell, ModuleSymmetry::Symmetry& symm);
    
    void get_irreps(UnitCell& ucell, ModuleSymmetry::Symmetry& symm);
};

} // namespace ModuleCell

#endif // QLIST_H
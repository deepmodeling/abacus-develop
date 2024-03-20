#ifndef CAL_R_OVERLAP_R_H
#define CAL_R_OVERLAP_R_H

#include <vector>
#include <map>
#include <set>

#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb21.h"

#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/vector3.h"
#include "module_base/abfs-vector3_order.h"
#include "module_base/ylm.h"
#include "single_R_io.h"


// output r_R matrix, added by Jingan
class cal_r_overlap_R
{

public:

    cal_r_overlap_R();
    ~cal_r_overlap_R();

    double kmesh_times = 4;
    double sparse_threshold = 1e-10;
    bool binary = false;

    void init(const Parallel_Orbitals &pv);
    void out_rR(const int &istep);
    void out_rR_other(const int &istep, const std::set<Abfs::Vector3_Order<int>> &output_R_coor);

private:
    void initialize_orb_table();
    void construct_orbs_and_orb_r();

    std::vector<int> iw2ia;
    std::vector<int> iw2iL;
    std::vector<int> iw2im;
    std::vector<int> iw2iN;
    std::vector<int> iw2it;

    ORB_table_phi MOT;
    ORB_gaunt_table MGT;

    Numerical_Orbital_Lm orb_r;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;

    std::map<size_t,
        std::map<size_t,
            std::map<size_t,
                std::map<size_t,
                    std::map<size_t,
                        std::map<size_t,
                        Center2_Orb::Orb11>>>>>> center2_orb11;
    
    std::map<size_t,
        std::map<size_t,
            std::map<size_t,
                std::map<size_t,
                    std::map<size_t,
                        std::map<size_t,
                        Center2_Orb::Orb21>>>>>> center2_orb21_r;

    const Parallel_Orbitals* ParaV;
    
};
#endif


#pragma once

#include <memory>
#include <vector>
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_cell/atom_spec.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint_type.h"
#include "big_grid.h"
#include "gint_atom.h"
#include "unitcell_info.h"
#include "localcell_info.h"
#include "divide_info.h"

#ifdef __CUDA
#include "batch_biggrid.h"
#include "source_lcao/module_gint/temp_gint/kernel/gint_gpu_vars.h"
#endif

namespace ModuleGint
{

class GintInfo
{
    public:
    // constructor
    GintInfo(
        int nbx, int nby, int nbz,
        int nmx, int nmy, int nmz,
        int startidx_bx, int startidx_by, int startidx_bz,
        int nbx_local, int nby_local, int nbz_local,
        const Numerical_Orbital* Phi,
        const UnitCell& ucell, Grid_Driver& gd);

    // getter functions
    const std::vector<std::shared_ptr<BigGrid>>& get_biggrids() { return biggrids_; }
    const std::vector<int>& get_trace_lo() const{ return trace_lo_; }
    int get_lgd() const { return lgd_; }
    int get_nat() const { return ucell_->nat; }        // return the number of atoms in the unitcell
    const UnitCell* get_ucell() const { return ucell_; }
    int get_local_mgrid_num() const { return localcell_info_->get_mgrids_num(); }
    double get_mgrid_volume() const { return meshgrid_info_->get_volume(); }

    //=========================================
    // functions about hcontainer
    //=========================================
    template <typename T>
    HContainer<T> get_hr(int npol = 1) const;
    
    private:
    // initialize the atoms
    void init_atoms_(int ntype, const Atom* atoms, const Numerical_Orbital* Phi);

    // initialize trace_lo_ and lgd_
    void init_trace_lo_(const UnitCell& ucell, const int nspin);

    // initialize the ijr_info
    void init_ijr_info_(const UnitCell& ucell, Grid_Driver& gd);

    const UnitCell* ucell_;

    // the unitcell information
    std::shared_ptr<const UnitCellInfo> unitcell_info_;

    // the biggrid information
    std::shared_ptr<const BigGridInfo> biggrid_info_;

    // the meshgrid information
    std::shared_ptr<const MeshGridInfo> meshgrid_info_;

    // the divide information
    std::shared_ptr<const DivideInfo> divide_info_;

    // the localcell information
    std::shared_ptr<const LocalCellInfo> localcell_info_;

    // the big grids on this processor
    std::vector<std::shared_ptr<BigGrid>> biggrids_;

    // the total atoms in the unitcell(include extended unitcell) on this processor
    // atoms[iat][Vec3i] is the atom with index iat in the unitcell with index Vec3i
    // Note: Since GintAtom does not implement a default constructor,
    // the map should not be accessed using [], but rather using the at function
    std::vector<std::map<Vec3i, GintAtom>> atoms_;

    // if the iat-th(global index) atom is in this processor, return true
    std::vector<bool> is_atom_in_proc_;

    // format for storing atomic pair information in hcontainer, used for initializing hcontainer
    std::vector<int> ijr_info_;

    // map the global index of atomic orbitals to local index
    std::vector<int> trace_lo_;
    
    // store the information about Numerical orbitals
    std::vector<Numerical_Orbital> orbs_;

    // total num of atomic orbitals on this proc
    int lgd_ = 0;

    #ifdef __CUDA
    public:
    std::vector<std::shared_ptr<BatchBigGrid>>& get_bgrid_batches() { return bgrid_batches_; };
    std::shared_ptr<const GintGpuVars> get_gpu_vars() const { return gpu_vars_; };
    int get_dev_id() const { return gpu_vars_->dev_id_; };
    int get_streams_num() const { return streams_num_; };
    
    private:
    void init_bgrid_batches_(int batch_size);
    std::vector<std::shared_ptr<BatchBigGrid>> bgrid_batches_;
    std::shared_ptr<const GintGpuVars> gpu_vars_;
    // More streams can improve parallelism and may speed up grid integration, at the cost of higher GPU memory usage.
    int streams_num_;
    #endif
};

} // namespace ModuleGint

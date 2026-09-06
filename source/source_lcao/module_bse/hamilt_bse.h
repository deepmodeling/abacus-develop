// Purpose: simplify the hamilt_casida class, and explore to go beyond Tamm-Damcoff Approximation.
// For simplicity, only ELPA solver with MPI parallization is implementated.
// Thus, instead of iterative solver such as Davidson, here matrix is constructed directly.

#pragma once
#include "bse_util.h"
#include "hamilt_bse_solver.h"
#include "molecular_lri.h"

#include "source_base/parallel_2d.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/hamilt.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/ao_to_mo.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/ri_benchmark/operator_ri_hartree.h"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/module_lr/utils/lr_util.h"

namespace BSE
{
template <typename T>
class HamiltBSE
{
public:
    std::vector<T> BSE_A_local, BSE_B_local;
    std::vector<T> VA_local, WA_local;
    std::vector<T> VB_local, WB_local;
    int ndim = 0; // dimension of BSE matrix
    Parallel_2D pA;
    /// @brief constructor for BSE_Matrix
    HamiltBSE(const int& nspin,
              const int& naos,
              const std::vector<int>& nocc,
              const std::vector<int>& nvirt,
              const UnitCell& ucell_in,
              const std::vector<double>& orb_cutoff_in,
              const Grid_Driver& gd_in,
              const psi::Psi<T>& psi_in,
              const psi::Psi<T>& psi_glb_in,
              const ModuleBase::matrix& eig_gw_in,
              MolecularLRI<T>& mo_lri_in,
              std::weak_ptr<LR::PotHxcLR> pot_in,
              const K_Vectors& kv_in,
              const std::vector<Parallel_2D>& pX_in,
              const Parallel_2D& pc_in,
              const Parallel_Orbitals& pmat_in,
              const std::vector<std::string>& spin_types_in,
              const std::string& tda,
              const bool bse_ri_hartree_in,
              const bool bse_mem_save_in,
              const int bse_continue_in,
              const bool out_bse_ab_in,
              const std::string& out_dir_in,
              const std::string& readin_dir_in,
              const int my_rank_in,
              const int nproc_in,
              const std::string& ri_hartree_benchmark_in = "none");
    /// @note these four functions shouldn't be called when `bse_mem_save` is true
    void cal_V_for_A();
    void cal_W_for_A();
    void cal_V_for_B();
    void cal_W_for_B();
    
    /// @brief initialize BSE matrix
    /// @note if `bse_mem_save` is true, V and W matrix will be added directly
    void init_bse_matrix(const bool is_full, const int& st_index);

    void tda_solver(const int& st_index, const int& nstates, double* ene_out, T* X_out);
    void full_solver(const int& st_index, const int& nstates, double* ene_out, T* X_out, T* Y_out);

    void cal_V_by_grid(bool is_A);
    void grid_calculation(hamilt::HContainer<T>& VR) const;
    
    inline void write_AB_matrix(const std::string& file, const int& prec, const T* ptr, const int& size1, const int& size2)
    {
        std::ofstream ofs(file);
        if (!ofs.is_open()){
            throw std::runtime_error("Cannot open file " + file);
        }
        ofs << file << "(Ry, transpose) with threshold " << prec << std::endl;
        ofs << std::setprecision(prec) << std::scientific;
        LR_Util::write_value(ofs, ptr, size1, size2);
        ofs.close();
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish writing " + file);
    }

    inline void read_AB_matrix(const std::string& file, T* ptr, const int& size1, const int& size2)
    {
        std::ifstream ifs(file);
        if (!ifs.is_open()){
            throw std::runtime_error("Cannot open file " + file);
        }
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip the first line
        LR_Util::read_value(ifs, ptr, size1, size2);
        ifs.close();
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "finish reading " + file);
    }

private:
    const int nspin;
    const int naos;
    const std::vector<int> nocc;
    const std::vector<int> nvirt;
    const UnitCell& ucell;
    const std::vector<double>& orb_cutoff;
    const Grid_Driver& gd;
    const psi::Psi<T>& psi_ks;
    const psi::Psi<T>& psi_ks_glb;
    const ModuleBase::matrix& eig_gw;
    MolecularLRI<T>& mo_lri;

    std::weak_ptr<LR::PotHxcLR> pot;
    const K_Vectors& kv;
    int nk = 1;
    const std::vector<Parallel_2D>& pX; // for tda, also pY for full
    const Parallel_2D& pc;
    const Parallel_Orbitals& pmat;
    const std::vector<std::string>& spin_types; // singlet, triplet, and rpa, ipa(independent particle approx)
    const bool bse_ri_hartree;
    const bool bse_mem_save;
    const int bse_continue;
    const bool out_bse_ab;
    const std::string out_dir;
    const std::string readin_dir;
    const int my_rank;
    const int nproc;
    const std::string ri_hartree_benchmark;

    std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_trans = nullptr;
};
} // namespace BSE

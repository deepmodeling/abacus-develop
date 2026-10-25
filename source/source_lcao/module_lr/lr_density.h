#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_H
#include "source_psi/psi.h"
#include "source_lcao/module_lr/dm_trans/dm_diff.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include <iosfwd>
#include <string>
class K_Vectors;
namespace LR
{
    template<typename T>
    class LR_Density
    {
        const UnitCell& ucell_;
        const K_Vectors& kv_;
        const Grid_Driver& gd_;
        const psi::Psi<T>& psi_ks_;
        const std::vector<double>& orb_cutoff_;
        const Parallel_Grid& pgrid_;
        const int nspin_;
        const std::vector<int>& nocc_;
        const std::vector<int>& nvirt_;
        const int nao_;
        int nk_;
        const std::vector<Parallel_2D>& pX_;
        const Parallel_2D& pc_;
        const Parallel_Orbitals& pmat_;
        const bool openshell_;
        std::vector<std::string> spintype_;
        /// cached aliases, read once here instead of at every call site below
        const int out_chg_precision_;
        const std::string global_out_dir_;

        void dm_to_density(module_dm::DensityMatrix<double, double>& dm, double** density);
        void dm_to_density(module_dm::DensityMatrix<std::complex<double>, std::complex<double>>& dm, double** density);

    public:
        LR_Density(const UnitCell& ucell,
            const K_Vectors& kv,
            const Grid_Driver& gd,
            const psi::Psi<T>& psi_ks,
            const std::vector<double>& orb_cutoff,
            const Parallel_Grid& pgrid,
            const int& nspin,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const int& nao,
            const std::vector<Parallel_2D>& pX,
            const Parallel_2D& pc,
            const Parallel_Orbitals& pmat,
            const bool openshell, const int output_precision, const std::string& output_directory);

        /// @brief calculate the electron density from the density matrix in 2d-block distribution
        void cal_eh_density_single_state(const T* const X_istate, const int ispin, double** density);

        void write_density_single_state(const double* const* const density, const std::string& filepath,
                                        std::ofstream& ofs_running);

        void output_eh_density_all_states(const T* const X, const int ispin, const int nstate,
                                          std::ofstream& ofs_running);
    };

}
#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_H

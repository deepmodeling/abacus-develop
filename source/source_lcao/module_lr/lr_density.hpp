#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_HPP
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_HPP
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_psi/psi.h"
#include "source_lcao/module_lr/dm_trans/dm_diff.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_io/module_output/cube_io.h"
#include "lr_amp.h"
#include <iosfwd>
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
        const int nk_;
        const std::vector<Parallel_2D>& pX_;
        const Parallel_2D& pc_;
        const Parallel_Orbitals& pmat_;
        const bool openshell_;
        const std::vector<std::string> spintype_;
        /// cached aliases, read once here instead of at every call site below
        const int out_chg_precision_ = PARAM.inp.out_chg[1];
        const std::string global_out_dir_ = PARAM.globalv.global_out_dir;

        inline void dm_to_density(module_dm::DensityMatrix<double, double>& dm, double** density)
        {
            ModuleBase::TITLE("LR_Density", "dm_to_density");
            ModuleGint::cal_gint_rho(dm.get_dmr_vec(), 1, density, false);
        }
        inline void dm_to_density(module_dm::DensityMatrix<std::complex<double>, std::complex<double>>& dm, double** density)
        {
            ModuleBase::TITLE("LR_Density", "dm_to_density");
            auto dm_to_density_real = [&](const char& part) -> void
                {
                    module_dm::DensityMatrix<std::complex<double>, double> dm_real(&pmat_, 1, kv_.kvec_d, nk_);
                    LR_Util::initialize_DMR<std::complex<double>, double>(dm_real, pmat_, ucell_, gd_, orb_cutoff_);
                    LR_Util::get_DMR_real_imag_part(dm, dm_real, part);
                    ModuleGint::cal_gint_rho(dm_real.get_dmr_vec(), 1, density, false);  // add-on
                };
            dm_to_density_real('R');
            dm_to_density_real('I');
        }

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
            const bool openshell = false) :
            ucell_(ucell), kv_(kv), gd_(gd), psi_ks_(psi_ks),
            orb_cutoff_(orb_cutoff), pgrid_(pgrid), nspin_(nspin), nk_(kv.get_nks() / nspin),
            nocc_(nocc), nvirt_(nvirt), nao_(nao), pX_(pX), pc_(pc), pmat_(pmat),
            openshell_(openshell), spintype_(openshell ? std::vector<std::string>({ "up", "down" }) : std::vector<std::string>({ "singlet", "triplet" }))
        {
        }

        /// @brief calculate the electron density from the density matrix in 2d-block distribution
        void cal_eh_density_single_state(const T* const X_istate, const int ispin, double** density)
        {
            ModuleBase::TITLE("LR_Density", "cal_eh_density_single_state");
            ModuleBase::GlobalFunc::ZEROS(density[0], this->pgrid_.get_nrxx());
            // 1. calculate the density matrix in AO basis
            auto c_spin = LR_Util::get_psi_spin(psi_ks_, ispin,nk_);
#ifdef __MPI
            const std::vector<ct::Tensor> dm_diff_k =
                cal_dm_diff_pblas<T>(X_istate, pX_[ispin], c_spin, pc_, nao_, nocc_[ispin], nvirt_[ispin], pmat_);
#else
            const std::vector<ct::Tensor> dm_diff_k =
                cal_dm_diff_blas<T>(X_istate, c_spin, nao_, nocc_[ispin], nvirt_[ispin]);
#endif
            // 2. calculate DM(R)
            module_dm::DensityMatrix<T, T> dm_diff=
                LR_Util::build_dm_from_dmk<T, T>(dm_diff_k,
                    this->pmat_, this->nk_, this->kv_.kvec_d, this->ucell_, this->gd_, this->orb_cutoff_);
            // 3. calculate electron density from DM(R)
            this->dm_to_density(dm_diff, density);
        }

        void write_density_single_state(const double* const* const density, const std::string& filepath,
                                        std::ofstream& ofs_running)
        {
            const std::string data_desc = "LR-TDDFT electron-hole density difference (electrons/Bohr^3)";
            ModuleIO::write_vdata_palgrid(pgrid_, density[0], 0, 1, 0, filepath, 0.0, &ucell_,
                                        this->out_chg_precision_, 0, false, true, ofs_running, data_desc);
        }

        void output_eh_density_all_states(const T* const X, const int ispin, const int nstate,
                                          std::ofstream& ofs_running)
        {
            ModuleBase::TITLE("LR_Density", "cal_eh_density_all_states");
            const int up_size = this->pX_[0].get_local_size();
            const int down_channel = openshell_ ? 1 : ispin;
            const int down_size = this->pX_[down_channel].get_local_size();
            double** density;
            LR_Util::_allocate_2order_nested_ptr(density, 1, pgrid_.get_nrxx());
            for (int istate = 0;istate < nstate;++istate)
            {
                const int offset = LR::electron_hole_offset(istate, nk_, up_size, down_size, ispin, openshell_);
                const T* const amplitudes = X + offset;
                this->cal_eh_density_single_state(amplitudes, ispin, density);
                const std::string filepath = this->global_out_dir_ + "LR_e-h_density_" + spintype_[ispin] + "_" + std::to_string(istate + 1) + ".cube";
                this->write_density_single_state(density, filepath, ofs_running);
            }
            LR_Util::_deallocate_2order_nested_ptr(density, 1);
        }
    };

}
#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_DENSITY_HPP

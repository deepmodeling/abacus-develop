#pragma once
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_psi/psi.h"
#include "source_lcao/module_lr/Grad/dm_diff/dm_diff.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_io/module_output/cube_io.h"
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
        
        inline void dm_to_density(elecstate::DensityMatrix<double, double>& dm, double** density)
        {
            ModuleBase::TITLE("LR_Density", "dm_to_density");
            ModuleGint::cal_gint_rho(dm.get_DMR_vector(), 1, density, false);
        }
        inline void dm_to_density(elecstate::DensityMatrix<std::complex<double>, std::complex<double>>& dm, double** density)
        {
            ModuleBase::TITLE("LR_Density", "dm_to_density");
            auto dm_to_density_real = [&](const char& part) -> void
                {
                    elecstate::DensityMatrix<std::complex<double>, double> dm_real(&pmat_, 1, kv_.kvec_d, nk_);
                    LR_Util::initialize_DMR<std::complex<double>, double>(dm_real, pmat_, ucell_, gd_, orb_cutoff_);
                    LR_Util::get_DMR_real_imag_part(dm, dm_real, part);
                    ModuleGint::cal_gint_rho(dm_real.get_DMR_vector(), 1, density, false);  // add-on
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
            const std::vector<ct::Tensor> dm_diff_k =
                cal_dm_diff_pblas<T>(X_istate, pX_[ispin], c_spin, pc_, nao_, nocc_[ispin], nvirt_[ispin], pmat_);
            // 2. calculate DM(R)
            elecstate::DensityMatrix<T, T> dm_diff=
                LR_Util::build_dm_from_dmk<T, T>(dm_diff_k,
                    this->pmat_, this->nk_, this->kv_.kvec_d, this->ucell_, this->gd_, this->orb_cutoff_);
            // 3. calculate electron density from DM(R)
            this->dm_to_density(dm_diff, density);
        }

        void write_density_single_state(const double* const* const density, const std::string& filepath)
        {
            ModuleIO::write_vdata_palgrid(pgrid_, density[0], 0, 1, 0, filepath, 0.0, &ucell_, PARAM.inp.out_chg[1], 0, false, true);
        }

        void output_eh_density_all_states(const T* const X, const int ispin, const int nstate)
        {
            ModuleBase::TITLE("LR_Density", "cal_eh_density_all_states");
            const int offset_per_state = openshell_ ?
                this->nk_ * (this->pX_[0].get_local_size() + this->pX_[1].get_local_size())
                : this->nk_ * this->pX_[ispin].get_local_size();
            double** density;
            LR_Util::_allocate_2order_nested_ptr(density, 1, pgrid_.get_nrxx());
            for (int istate = 0;istate < nstate;++istate)
            {
                int offset = istate * offset_per_state;
                if (openshell_)
                    offset += ispin * this->pX_[0].get_local_size();
                this->cal_eh_density_single_state(X + offset, ispin, density);
                const std::string filepath = PARAM.globalv.global_out_dir + "LR_e-h_density_" + spintype_[ispin] + "_" + std::to_string(istate + 1) + ".cube";
                this->write_density_single_state(density, filepath);
            }
            LR_Util::_deallocate_2order_nested_ptr(density, 1);
        }
    };

}
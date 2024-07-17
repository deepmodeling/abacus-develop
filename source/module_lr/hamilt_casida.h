#pragma once
#include "module_hamilt_general/hamilt.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/operator_casida/operator_lr_diag.h"
#include "module_lr/operator_casida/operator_lr_hxc.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#ifdef __EXX
#include "module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{
    template<typename T>
    class HamiltCasidaLR : public hamilt::Hamilt<T, base_device::DEVICE_CPU>
    {
    public:
        template<typename TGint>
        HamiltCasidaLR(std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const int& nocc,
            const int& nvirt,
            const UnitCell& ucell_in,
            Grid_Driver& gd_in,
            const psi::Psi<T>* psi_ks_in,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            Exx_LRI<T>* exx_lri_in,
            const double& exx_alpha,
#endif 
            TGint* gint_in,
            std::shared_ptr<PotHxcLR> pot_in,
            const K_Vectors& kv_in,
            Parallel_2D* pX_in,
            Parallel_2D* pc_in,
            Parallel_Orbitals* pmat_in) : nocc(nocc), nvirt(nvirt), pX(pX_in), nk(kv_in.get_nks() / nspin)
        {
            ModuleBase::TITLE("HamiltCasidaLR", "HamiltCasidaLR");
            this->classname = "HamiltCasidaLR";
            this->DM_trans.resize(1);
            this->DM_trans[0] = new elecstate::DensityMatrix<T, T>(&kv_in, pmat_in, nspin);
            // add the diag operator  (the first one)
            this->ops = new OperatorLRDiag<T>(eig_ks, pX_in, nk, nocc, nvirt);
            //add Hxc operator
            OperatorLRHxc<T>* lr_hxc = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks_in,
                this->DM_trans, gint_in, &(*pot_in), ucell_in, gd_in, kv_in, pX_in, pc_in, pmat_in);
            this->ops->add(lr_hxc);
#ifdef __EXX
            if (xc_kernel == "hf" || xc_kernel == "hse")
            {   //add Exx operator
                hamilt::Operator<T>* lr_exx = new OperatorLREXX<T>(nspin, naos, nocc, nvirt, ucell_in, psi_ks_in,
                    this->DM_trans, exx_lri_in, kv_in, pX_in, pc_in, pmat_in, exx_alpha);
                this->ops->add(lr_exx);
            }
#endif
        }
        ~HamiltCasidaLR()
        {
            if (this->ops != nullptr)
            {
                delete this->ops;
            }
            for (auto& d : this->DM_trans)delete d;
        };

        hamilt::HContainer<T>* getHR() { return this->hR; }

        virtual std::vector<T> matrix() override;

    private:
        int nocc;
        int nvirt;
        int nk;
        Parallel_2D* pX = nullptr;
        T one();
        hamilt::HContainer<T>* hR = nullptr;
        /// transition density matrix in AO representation
        /// Hxc only: size=1, calculate on the same address for each bands
        /// Hxc+Exx: size=nbands, store the result of each bands for common use
        std::vector<elecstate::DensityMatrix<T, T>*> DM_trans;
    };
}
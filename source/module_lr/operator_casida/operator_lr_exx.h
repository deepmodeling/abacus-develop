#ifdef __EXX
#pragma once
#include "module_hamilt_general/operator.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_ri/Exx_LRI.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{

    /// @brief  Hxc part of A operator
    template<typename T = double>
    class OperatorLREXX : public hamilt::Operator<T, base_device::DEVICE_CPU>
    {
        using TA = int;
        static const size_t Ndim = 3;
        using TC = std::array<int, Ndim>;
        using TAC = std::pair<TA, TC>;

    public:
        OperatorLREXX(const int& nspin,
            const int& naos,
            const int& nocc,
            const int& nvirt,
            const UnitCell& ucell_in,
            const psi::Psi<T>& psi_ks_in,
            std::unique_ptr<elecstate::DensityMatrix<T, T>>& DM_trans_in,
            // HContainer<double>* hR_in,
            std::weak_ptr<Exx_LRI<T>> exx_lri_in,
            const K_Vectors& kv_in,
            const Parallel_2D& pX_in,
            const Parallel_2D& pc_in,
            const Parallel_Orbitals& pmat_in,
            const double& alpha = 1.0,
            const std::vector<int>& aims_nbasis = {})
            : nspin(nspin), naos(naos), nocc(nocc), nvirt(nvirt), nk(kv_in.get_nks() / nspin),
            psi_ks(psi_ks_in), DM_trans(DM_trans_in), exx_lri(exx_lri_in), kv(kv_in),
            pX(pX_in), pc(pc_in), pmat(pmat_in), ucell(ucell_in), alpha(alpha),
            aims_nbasis(aims_nbasis)
        {
            ModuleBase::TITLE("OperatorLREXX", "OperatorLREXX");
            this->cal_type = hamilt::calculation_type::lcao_exx;
            this->is_first_node = false;

            // reduce psi_ks for later use
            this->psi_ks_full.resize(this->nk, nocc + nvirt, this->naos);
            for (int ik = 0;ik < nk;++ik)
            {
                LR_Util::gather_2d_to_full(this->pc, &this->psi_ks(ik, 0, 0), &this->psi_ks_full(ik, 0, 0), false, this->naos, nocc + nvirt);
            }

            // get cells in BvK supercell
            const TC period = RI_Util::get_Born_vonKarmen_period(kv_in);
            this->BvK_cells = RI_Util::get_Born_von_Karmen_cells(period);

            this->allocate_Ds_onebase();
            this->exx_lri.lock()->Hexxs.resize(1);
        };

        void init(const int ik_in) override {};

        virtual void act(const int nbands,
                         const int nbasis,
                         const int npol,
                         const T* psi_in,
                         T* hpsi,
                         const int ngk_ik = 0,
                         const bool is_first_node = false) const override;

    private:
        //global sizes
        const int nspin = 1;
        const int naos = 1;
        const int nocc = 1;
        const int nvirt = 1;
        const int nk = 1;  ///< number of k-points
        const double alpha = 1.0;   //(allow non-ref constant)
        const bool cal_dm_trans = false;
        const bool tdm_sym = false; ///< whether transition density matrix is symmetric
        const K_Vectors& kv;
        /// ground state wavefunction
        const psi::Psi<T>& psi_ks = nullptr;
        psi::Psi<T> psi_ks_full;
        const std::vector<int> aims_nbasis={};    ///< number of basis functions for each type of atom in FHI-aims

        /// transition density matrix 
        std::unique_ptr<elecstate::DensityMatrix<T, T>>& DM_trans;

        /// density matrix of a certain (i, a, k), with full naos*naos size for each key
        /// D^{iak}_{\mu\nu}(k): 1/N_k * c^*_{ak,\mu} c_{ik,\nu}
        /// D^{iak}_{\mu\nu}(R): D^{iak}_{\mu\nu}(k)e^{-ikR}
        // elecstate::DensityMatrix<T, double>* DM_onebase;
        mutable std::map<TA, std::map<TAC, RI::Tensor<T>>> Ds_onebase;

        // cells in the Born von Karmen supercell (direct)
        std::vector<std::array<int, Ndim>> BvK_cells;

        /// transition hamiltonian in AO representation
        // hamilt::HContainer<double>* hR = nullptr;

        /// C, V tensors of RI, and LibRI interfaces
        /// gamma_only: T=double, Tpara of exx (equal to Tpara of Ds(R) ) is also double 
        ///.multi-k: T=complex<double>, Tpara of exx here must be complex, because Ds_onebase is complex
        /// so TR in DensityMatrix and Tdata in Exx_LRI are all equal to T
        std::weak_ptr<Exx_LRI<T>> exx_lri;

        const UnitCell& ucell;

        ///parallel info
        const Parallel_2D& pc;
        const Parallel_2D& pX;
        const Parallel_Orbitals& pmat;


        // allocate Ds_onebase
        void allocate_Ds_onebase();

        void cal_DM_onebase(const int io, const int iv, const int ik) const;

    };
}
#endif
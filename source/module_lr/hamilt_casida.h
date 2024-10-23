#pragma once
#include <typeinfo>
#include "module_hamilt_general/hamilt.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/operator_casida/operator_lr_diag.h"
#include "module_lr/operator_casida/operator_lr_hxc.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_lr/dm_trans/dm_trans.h"
#ifdef __EXX
#include "module_lr/operator_casida/operator_lr_exx.h"
#include "module_lr/ri_benchmark/operator_ri_hartree.h"
#include "module_ri/LRI_CV_Tools.h"
#endif
namespace LR
{
    template<typename T>
    class HamiltLR
    {
    public:
        template<typename TGint>
        HamiltLR(std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const UnitCell& ucell_in,
            const std::vector<double>& orb_cutoff,
            Grid_Driver& gd_in,
            const psi::Psi<T>& psi_ks_in,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            std::weak_ptr<Exx_LRI<T>> exx_lri_in,
            const double& exx_alpha,
#endif 
            TGint* gint_in,
            std::weak_ptr<PotHxcLR> pot_in,
            const K_Vectors& kv_in,
            const std::vector<Parallel_2D>& pX_in,
            const Parallel_2D& pc_in,
            const Parallel_Orbitals& pmat_in,
            const std::string& spin_type,
            const std::string& ri_hartree_benchmark = "none",
            const std::vector<int>& aims_nbasis = {}) : nspin(nspin), nocc(nocc), nvirt(nvirt), pX(pX_in), nk(kv_in.get_nks() / nspin)
        {
            ModuleBase::TITLE("HamiltLR", "HamiltLR");
            if (ri_hartree_benchmark != "aims") { assert(aims_nbasis.empty()); }
            // always use nspin=1 for transition density matrix
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat_in, 1, kv_in.kvec_d, nk);
            if (ri_hartree_benchmark == "none") { LR_Util::initialize_DMR(*this->DM_trans, pmat_in, ucell_in, gd_in, orb_cutoff); }
            // this->DM_trans->init_DMR(&gd_in, &ucell_in); // too large due to not restricted by orb_cutoff

            // add the diag operator  (the first one)
            this->ops = new OperatorLRDiag<T>(eig_ks.c, pX[0], nk, nocc[0], nvirt[0]);
            //add Hxc operator
#ifdef __EXX
            using TAC = std::pair<int, std::array<int, 3>>;
            using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;
            const std::string& dir = PARAM.globalv.global_readin_dir;
            TLRI Cs_read; 
            TLRI Vs_read; 
#ifdef __DEBUG
            // TLRI Vs_compare = LRI_CV_Tools::read_Vs_abf<T>(dir + "Vs");
            // LRI_CV_Tools::write_Vs_abf(Vs_read, "Vs_read_from_coulomb");
            // LRI_CV_Tools::write_Cs_ao(Cs_read, "Cs_ao_read"); // ensure Cs_ao is read correctly
            // assert(RI_Benchmark::compare_Vs(Vs_read, Vs_compare));
#endif
            if (ri_hartree_benchmark != "none")
            {
#ifdef __EXX
                if (spin_type == "singlet")
                {
                    if (ri_hartree_benchmark == "aims") 
                    { 
                        Cs_read = LRI_CV_Tools::read_Cs_ao<T>(dir + "Cs_data_0.txt");
                        Vs_read = RI_Benchmark::read_coulomb_mat_general<T>(dir + "coulomb_mat_0.txt", Cs_read); 
                    }
                    else if (ri_hartree_benchmark == "abacus")
                    {
                        Cs_read = LRI_CV_Tools::read_Cs_ao<T>(dir + "Cs");
                        Vs_read = LRI_CV_Tools::read_Vs_abf<T>(dir + "Vs");
                    }
                    if (!std::set<std::string>({ "rpa", "hf" }).count(xc_kernel)) { throw std::runtime_error("ri_hartree_benchmark is only supported for xc_kernel rpa and hf"); }
                    RI_Benchmark::OperatorRIHartree<T>* ri_hartree_op
                        = new RI_Benchmark::OperatorRIHartree<T>(ucell_in, naos, nocc[0], nvirt[0], psi_ks_in,
                            Cs_read, Vs_read, ri_hartree_benchmark == "aims", aims_nbasis);
                    this->ops->add(ri_hartree_op);
                }
                else if (spin_type == "triplet") { std::cout << "f_Hxc based on grid integral is not needed." << std::endl; }
#else
                ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
            }
            else
#endif
            {
                OperatorLRHxc<T>* lr_hxc = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks_in,
                    this->DM_trans, gint_in, pot_in, ucell_in, orb_cutoff, gd_in, kv_in, pX_in, pc_in, pmat_in);
                this->ops->add(lr_hxc);
            }
#ifdef __EXX
            if (xc_kernel == "hf" || xc_kernel == "hse")
            {   //add Exx operator
                if (ri_hartree_benchmark != "none" && spin_type == "singlet")
                {
                    exx_lri_in.lock()->reset_Cs(Cs_read);
                    exx_lri_in.lock()->reset_Vs(Vs_read);
                }
                // std::cout << "exx_alpha=" << exx_alpha << std::endl; // the default value of exx_alpha is 0.25 when dft_functional is pbe or hse
                hamilt::Operator<T>* lr_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell_in, psi_ks_in,
                    this->DM_trans, exx_lri_in, kv_in, pX_in[0], pc_in, pmat_in,
                    xc_kernel == "hf" ? 1.0 : exx_alpha, //alpha
                    aims_nbasis);
                this->ops->add(lr_exx);
            }
#endif

            this->cal_dm_trans = [&, this](const int& is, const T* X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks_in, is, nk);
#ifdef __MPI
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X, pX[is], psi_ks_is, pc_in, naos, nocc[is], nvirt[is], pmat_in);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat_in);
#else
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is]);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                    // LR_Util::print_tensor<T>(dm_trans_2d[0], "dm_trans_2d[0]", &pmat_in);
                    // tensor to vector, then set DMK
                    for (int ik = 0;ik < nk;++ik) { this->DM_trans->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                };
        }
        ~HamiltLR() { delete this->ops; }

        std::vector<T> matrix()const;

        void hPsi(const T* psi_in, T* hpsi, const int ld_psi, const int& nband) const
        {
            assert(ld_psi == nk * pX[0].get_local_size());
            for (int ib = 0;ib < nband;++ib)
            {
                const int offset = ib * ld_psi;
                this->cal_dm_trans(0, psi_in + offset);  // calculate transition density matrix here
                hamilt::Operator<T>* node(this->ops);
                while (node != nullptr)
                {
                    node->act(/*nband=*/1, ld_psi, /*npol=*/1, psi_in + offset, hpsi + offset);
                    node = (hamilt::Operator<T>*)(node->next_op);
                }
            }
        }

        void global2local(T* lvec, const T* gvec, const int& nband) const
        {
            const int npairs = nocc[0] * nvirt[0];
            for (int ib = 0;ib < nband;++ib)
            {
                const int loffset_b = ib * nk * pX[0].get_local_size();
                const int goffset_b = ib * nk * npairs;
                for (int ik = 0;ik < nk;++ik)
                {
                    const int loffset = loffset_b + ik * pX[0].get_local_size();
                    const int goffset = goffset_b + ik * npairs;
                    for (int lo = 0;lo < pX[0].get_col_size();++lo)
                    {
                        const int go = pX[0].local2global_col(lo);
                        for (int lv = 0;lv < pX[0].get_row_size();++lv)
                        {
                            const int gv = pX[0].local2global_row(lv);
                            lvec[loffset + lo * pX[0].get_row_size() + lv] = gvec[goffset + go * nvirt[0] + gv];
                        }
                    }
                }
            }
        }

    private:
        const std::vector<int>& nocc;
        const std::vector<int>& nvirt;
        const int nspin = 1;
        const int nk = 1;
        const bool tdm_sym = false;     ///< whether to symmetrize the transition density matrix
        const std::vector<Parallel_2D>& pX;
        T one()const;
        /// transition density matrix in AO representation
        /// calculate on the same address for each bands, and commonly used by all the operators
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_trans;

        /// first node operator, add operations from each operators
        hamilt::Operator<T, base_device::DEVICE_CPU>* ops = nullptr;

        std::function<void(const int&, const T*)> cal_dm_trans;
    };
}

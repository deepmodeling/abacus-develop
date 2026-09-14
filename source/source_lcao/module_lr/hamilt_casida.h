#pragma once
#include <typeinfo>
#include "source_hamilt/hamilt.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_diag.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_hxc.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#include "source_lcao/module_lr/ri_benchmark/operator_ri_hartree.h"
#include "source_lcao/module_ri/lri_cv_tools.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#endif
namespace LR
{
    template<typename T>
    class HamiltLR
    {
    public:
      HamiltLR(std::string& xc_kernel,
               const int& nspin,
               const int& naos,
               const std::vector<int>& nocc,
               const std::vector<int>& nvirt,
               const UnitCell& ucell_in,
               const std::vector<double>& orb_cutoff,
               const Grid_Driver& gd_in,
               const psi::Psi<T>& psi_ks_in,
               const ModuleBase::matrix& eig_ks,
#ifdef __EXX
               std::weak_ptr<Exx_LRI<T>> exx_lri_in,
               const double& exx_alpha,
#endif
               std::weak_ptr<PotHxcLR> pot_in,
               const K_Vectors& kv_in,
               const std::vector<Parallel_2D>& pX_in,
               const Parallel_2D& pc_in,
               const Parallel_Orbitals& pmat_in,
               const std::string& spin_type,
               const std::string& in_dir,
               const std::string& out_dir,
               const std::string& ri_hartree_benchmark = "none",
               const std::vector<int>& aims_nbasis = {})
          : nspin(nspin), nocc(nocc), nvirt(nvirt), pX(pX_in), nk(kv_in.get_nks() / nspin)
        {
            ModuleBase::TITLE("HamiltLR", "HamiltLR");
            if (ri_hartree_benchmark != "aims" && ri_hartree_benchmark !="aims-librpa") { assert(aims_nbasis.empty()); }
            // always use nspin=1 for transition density matrix
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat_in, 1, kv_in.kvec_d, nk);
            if (ri_hartree_benchmark == "none") { LR_Util::initialize_DMR(*this->DM_trans, pmat_in, ucell_in, gd_in, orb_cutoff); }
            // this->DM_trans->init_DMR(&gd_in, &ucell_in); // too large due to not restricted by orb_cutoff

            // 1.add the diag operator  (the first one)
            this->ops = new OperatorLRDiag<T>(eig_ks.c, pX[0], nk, nocc[0], nvirt[0]);
            // 2.add Hxc operator
#ifdef __EXX
            using TAC = std::pair<int, std::array<int, 3>>;
            using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;
            TLRI Cs_read; 
            TLRI Vs_read; 
#ifdef __DEBUG
            // TLRI Vs_compare = LRI_CV_Tools::read_Vs_abf<T>(in_dir + "Vs");
            // LRI_CV_Tools::write_Vs_abf(Vs_read, "Vs_read_from_coulomb");
            // LRI_CV_Tools::write_Cs_ao(Cs_read, "Cs_ao_read"); // ensure Cs_ao is read correctly
            // assert(RI_Benchmark::compare_Vs(Vs_read, Vs_compare));
#endif
            if (ri_hartree_benchmark != "none")
            {
#ifdef __EXX
                if (spin_type == "singlet")
                {
                    int use_fine_kgrid = 0;
                    if (ri_hartree_benchmark == "aims" || ri_hartree_benchmark == "aims-librpa")
                    {
                        LR_IO::RI_kRlist kRlist (ucell_in,
                                                const_cast<K_Vectors*>(&kv_in),
                                                nspin, in_dir, out_dir, use_fine_kgrid);
                        // though C and V are real, here still use <T> to multiply with psi
                        Cs_read = LRI_CV_Tools::read_Cs_ao_all<T>(in_dir);
                        Vs_read = LR_IO::read_coulomb_mat_general_k<T,T>(in_dir, Cs_read, kRlist);
                    }
                    else if (ri_hartree_benchmark == "abacus")
                    {
                        Cs_read = LRI_CV_Tools::read_Cs_ao<T>(in_dir + "Cs");
                        Vs_read = LRI_CV_Tools::read_Vs_abf<T>(in_dir + "Vs");
                    }
                    else if (ri_hartree_benchmark == "abacus-librpa")
                    {
                        LR_IO::RI_kRlist kRlist (ucell_in,
                                                const_cast<K_Vectors*>(&kv_in),
                                                nspin, in_dir, out_dir, use_fine_kgrid);
                        Cs_read = LRI_CV_Tools::read_Cs_ao_all<T>(in_dir);
                        Vs_read = LR_IO::read_coulomb_mat_k<T,T>(in_dir, Cs_read, kRlist);
                    }
                    if (!std::set<std::string>({ "rpa", "hf"}).count(xc_kernel)) {
                        throw std::runtime_error("ri_hartree_benchmark is only supported for xc_kernel = rpa, hf");
                    }
                    RI_Benchmark::OperatorRIHartree<T>* ri_hartree_op
                        = new RI_Benchmark::OperatorRIHartree<T>(ucell_in, naos, nocc[0], nvirt[0], psi_ks_in,
                            Cs_read, Vs_read);
                    this->ops->add(ri_hartree_op);
                }
                else if (spin_type == "triplet") { std::cout << "Hatree term is not needed for S2:triplet." << std::endl; }
#else
                ModuleBase::WARNING_QUIT("ESolver_LR", "RI benchmark is only supported when compile with LibRI.");
#endif
            }
            else
#endif
            {
                OperatorLRHxc<T>* lr_hxc = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks_in,
                    this->DM_trans, pot_in, ucell_in, orb_cutoff, gd_in, kv_in, pX_in, pc_in, pmat_in);
                this->ops->add(lr_hxc);
            }
#ifdef __EXX// 3.add Exx operator
            if (xc_kernel == "hf" || xc_kernel == "hse")
            {
                if (ri_hartree_benchmark != "none" && spin_type == "singlet")
                {
                    exx_lri_in.lock()->reset_Cs(Cs_read);
                    exx_lri_in.lock()->reset_Vs(Vs_read);
                }
                // std::cout << "exx_alpha=" << exx_alpha << std::endl; // the default value of exx_alpha is 0.25 when dft_functional is pbe or hse
                hamilt::Operator<T>* lr_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell_in, psi_ks_in,
                    this->DM_trans, exx_lri_in, kv_in, pX_in[0], pc_in, pmat_in,
                    (xc_kernel == "hf") ? 1.0 : exx_alpha);
                this->ops->add(lr_exx);
            }
#endif

            this->cal_dm_trans = [&, this](const int& is, const T* const X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks_in, is, nk);
#ifdef __MPI
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X, pX[is], psi_ks_is, pc_in, naos, nocc[is], nvirt[is], pmat_in, (T)1.0 / (T)nk);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat_in);
#else
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is], (T)1.0 / (T)nk);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                    // LR_Util::print_tensor<T>(dm_trans_2d[0], "dm_trans_2d[0]", &pmat_in);
                    // tensor to vector, then set DMK
                    for (int ik = 0;ik < nk;++ik) { this->DM_trans->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                };
        }
        ~HamiltLR() { delete this->ops; }

        std::vector<T> matrix()const;

        void hPsi(const T* const psi_in, T* const hpsi, const int ld_psi, const int& nband) const
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

        // void global2local(T* lvec, const T* gvec, const int& nband) const
        // {
        //     const int npairs = nocc[0] * nvirt[0];
        //     for (int ib = 0;ib < nband;++ib)
        //     {
        //         const int loffset_b = ib * nk * pX[0].get_local_size();
        //         const int goffset_b = ib * nk * npairs;
        //         for (int ik = 0;ik < nk;++ik)
        //         {
        //             const int loffset = loffset_b + ik * pX[0].get_local_size();
        //             const int goffset = goffset_b + ik * npairs;
        //             for (int lo = 0;lo < pX[0].get_col_size();++lo)
        //             {
        //                 const int go = pX[0].local2global_col(lo);
        //                 for (int lv = 0;lv < pX[0].get_row_size();++lv)
        //                 {
        //                     const int gv = pX[0].local2global_row(lv);
        //                     lvec[loffset + lo * pX[0].get_row_size() + lv] = gvec[goffset + go * nvirt[0] + gv];
        //                 }
        //             }
        //         }
        //     }
        // }

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

        std::function<void(const int&, const T* const)> cal_dm_trans;
    };
}

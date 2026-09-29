#ifdef __EXX
#include "operator_lr_exx.h"
#include <cstdlib>
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_lcao/module_lr/ri_benchmark/ri_benchmark.h"
#include "source_lcao/module_lr/dm_band/dm_band.h"
namespace LR
{
    template<typename T>
    void OperatorLREXX<T>::allocate_Ds_onebase()
    {
        ModuleBase::TITLE("OperatorLREXX", "allocate_Ds_onebase");
        for (int iat1 = 0;iat1 < ucell.nat;++iat1) {
            const int it1 = ucell.iat2it[iat1];
            for (int iat2 = 0;iat2 < ucell.nat;++iat2) {
                const int it2 = ucell.iat2it[iat2];
                for (auto cell : this->BvK_cells) {
                    this->Ds_onebase[iat1][std::make_pair(iat2, cell)] = 
                        RI::Tensor<T>({ static_cast<size_t>(ucell.atoms[it1].nw),  static_cast<size_t>(ucell.atoms[it2].nw) });
                }
            }
        }
    }

    template<>
    void OperatorLREXX<double>::cal_DM_onebase(const int io, const int iv, const int ik) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
        switch (this->dm_pq_)
        {
        case MO_TO_AO_TYPE::CC_vo:
        {
            // Co Cv
            DMBand<double>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->psi_ks_full)
                .cal_dm_band(io, nocc + iv, ik, this->Ds_onebase, 1.0, this->aims_nbasis, this->aims_nbasis);
            break;
        }
        case MO_TO_AO_TYPE::CXC:
        {
            // term1: Co -> CvX, i.e. [CvX] Cv
            DMBand<double> dm_band1(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->cvx_full, this->psi_ks_full);
            dm_band1.eval(io, nocc + iv, ik);
            // term2: Cv -> CoX^T, i.e. Co [CoX^T]
            DMBand<double> dm_band2(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->coxt_full);
            dm_band2.eval(io, iv, ik);
            this->Ds_onebase = (dm_band1 - dm_band2).get_data();
            break;
        }
        case MO_TO_AO_TYPE::CC_oo:
        {
            // Cv -> CoX^T, i.e. C_o [C_oX^T]
            DMBand<double>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->psi_ks_full)
                .cal_dm_band(io, iv, ik, this->Ds_onebase, 1.0);
            break;
        }
        case MO_TO_AO_TYPE::CXC_o:
        {
            // Cv -> CoX^T, i.e. C_o [C_oX^T] (the same as CXC term2 but with positive sign)
            DMBand<double>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->coxt_full)
                .cal_dm_band(io, iv, ik, this->Ds_onebase);
            break;
        }
        default:
            break;
        }
    }

    template<>
    void OperatorLREXX<std::complex<double>>::cal_DM_onebase(const int io, const int iv, const int ik) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
        switch (this->dm_pq_)
        {
        case MO_TO_AO_TYPE::CC_vo:
        {
            // Cv Co^*
            DMBand<std::complex<double>>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->psi_ks_full)
                .cal_dm_band(nocc + iv, io, ik, this->Ds_onebase, 1.0, this->aims_nbasis, this->aims_nbasis);
            break;
        }
        case MO_TO_AO_TYPE::CXC:
        {
            // term1: Co -> CvX, i.e. Cv [CvX]^*
            DMBand<std::complex<double>> dm_band1(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->cvx_full);
            dm_band1.eval(nocc + iv, io, ik);
            // term2: Cv -> CoX^T, i.e. [CoX^T] Co^*
            DMBand<std::complex<double>> dm_band2(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->coxt_full, this->psi_ks_full);
            dm_band2.eval(iv, io, ik);
            this->Ds_onebase = (dm_band1 - dm_band2).get_data();
            break;
        }
        case MO_TO_AO_TYPE::CC_oo:
        {
            // Co Co^*
            // ! note: iv traverses nocc here, i.e. io=i, iv=j
            DMBand<std::complex<double>>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->psi_ks_full, this->psi_ks_full)
                .cal_dm_band(io, iv, ik, this->Ds_onebase, 1.0, this->aims_nbasis, this->aims_nbasis);
            break;
        }
        case MO_TO_AO_TYPE::CXC_o:
        {
            // Cv -> CoX^T, i.e. [C_oX^T]  C_o^* (the same as CXC term2 but with positive sign)
            DMBand<std::complex<double>>(ucell, pmat, this->kv.kvec_c, this->BvK_cells, this->coxt_full, this->psi_ks_full)
                .cal_dm_band(iv, io, ik, this->Ds_onebase);
            break;
        }
        default:
            break;
        }
    }

    template<typename T>
    void OperatorLREXX<T>::act(const int nbands, 
                               const int nbasis, 
                               const int npol, 
                               const T* psi_in, 
                               T* hpsi, 
                               const int ngk_ik, 
                               const bool is_first_node)const
    {
        ModuleBase::TITLE("OperatorLREXX", "act");
        ModuleBase::timer::start("OperatorLREXX", "act");

        // convert parallel info to LibRI interfaces
        std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(ucell,this->pmat);

        // suppose Cs，Vs, have already been calculated in the ion-step of ground state
        // and DM_trans has been calculated in hPsi() outside.

        // 1. set_Ds (once)
        // convert to vector<T*> for the interface of RI_2D_Comm::split_m2D_ktoR (interface will be unified to ct::Tensor)
        std::vector<std::vector<T>> DMk_trans_vector = this->DM_trans->get_dmk_vec();
        // assert(DMk_trans_vector.size() == nk);
        std::vector<const std::vector<T>*> DMk_trans_pointer(nk);
        for (int ik = 0;ik < nk;++ik) { DMk_trans_pointer[ik] = &DMk_trans_vector[ik]; }
        // if multi-k, DM_trans(TR=double) -> Ds_trans(TR=T=complex<double>)
        std::vector<std::map<TA, std::map<TAC, RI::Tensor<T>>>> Ds_trans =
            // aims_nbasis.empty() ? // ucell.nw is updated, abandoned 25-05-23
            RI_2D_Comm::split_m2D_ktoR<T>(ucell,this->kv, DMk_trans_pointer, this->pmat, 1);
            //: RI_Benchmark::split_Ds(DMk_trans_vector, aims_nbasis, ucell); //0.5 will be multiplied
        // LR_Util::print_CV(Ds_trans[0], "Ds_trans in OperatorLREXX", 1e-10);
        
        // 2. cal_Hs
        ModuleBase::timer::start("OperatorLREXX", "cal_Hs");
        auto lri = this->exx_lri.lock();
        lri->exx_lri.set_Ds(std::move(Ds_trans[0]), lri->info.dm_threshold);
        lri->exx_lri.cal_Hs();
        lri->Hexxs[0] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
            lri->mpi_comm, std::move(lri->exx_lri.Hs), std::get<0>(judge[0]), std::get<1>(judge[0]));
        lri->post_process_Hexx(lri->Hexxs[0]);
        // LR_Util::print_CV(lri->Hexxs[0], "Hexxs in OperatorLREXX", 1e-10);
        ModuleBase::timer::end("OperatorLREXX", "cal_Hs");
        // LR_Util::print_CV(lri->Hexxs[0], "Hexx in OperatorLREXX after post_process_Hexx", 1e-10);

        // 3. set [AX]_iak = DM_onbase * Hexxs for each occ-virt pair and each k-point
        // caution: parrallel
        ModuleBase::timer::start("OperatorLREXX", "cal_energy");
        if (this->dm_pq_ == MO_TO_AO_TYPE::CXC || this->dm_pq_ == MO_TO_AO_TYPE::CXC_o)
        {
            this->cal_coxt_cvx(psi_in);
        }
        const int nrow_global = (this->dm_pq_ == MO_TO_AO_TYPE::CC_oo) ? this->nocc : this->nvirt;
        for (int io = 0;io < this->nocc;++io)
        {
            for (int iv = 0;iv < nrow_global;++iv)
            {
                for (int ik = 0;ik < nk;++ik)
                {
                    const int xstart_bk = ik * pX.get_local_size();
                    this->cal_DM_onebase(io, iv, ik);       //set Ds_onebase for all e-h pairs (not only on this processor)
                    // LR_Util::print_CV(Ds_onebase, "Ds_onebase of occ " + std::to_string(io) + ", virtual " + std::to_string(iv) + " in OperatorLREXX", 1e-10);
                    const T& ene = 2 * alpha * //minus for exchange and Hartree-to-Ry are already considered in `post_process_Hexx`, 2 for canceling 0.5 in split_m2D_ktoR(nspin=1)
                        lri->exx_lri.post_2D.cal_energy(this->Ds_onebase, lri->Hexxs[0]);
                    if (this->pX.in_this_processor(iv, io))
                    {
                        hpsi[xstart_bk + this->pX.global2local_col(io) * this->pX.get_row_size() + this->pX.global2local_row(iv)] += ene;
                    }
                }
            }
        }
        ModuleBase::timer::end("OperatorLREXX", "cal_energy");
        ModuleBase::timer::end("OperatorLREXX", "act");
    }
    template class OperatorLREXX<double>;
    template class OperatorLREXX<std::complex<double>>;
}
#endif
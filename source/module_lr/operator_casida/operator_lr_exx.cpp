#ifdef __EXX
#include "operator_lr_exx.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_print.h"
#include "module_lr/ri_benchmark/ri_benchmark.h"
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
                    this->Ds_onebase[iat1][std::make_pair(iat2, cell)] = aims_nbasis.empty() ?
                        RI::Tensor<T>({ static_cast<size_t>(ucell.atoms[it1].nw),  static_cast<size_t>(ucell.atoms[it2].nw) }) :
                        RI::Tensor<T>({ static_cast<size_t>(aims_nbasis[it1]),  static_cast<size_t>(aims_nbasis[it2]) });
                }
            }
        }
    }

    template<>
    void OperatorLREXX<double>::cal_DM_onebase(const int io, const int iv, const int ik) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
        // NOTICE: DM_onebase will be passed into `cal_energy` interface and conjugated by "zdotc". 
        // So the formula should be the same as RHS. instead of LHS of the A-matrix, 
        // i.e. c1v · conj(c2o) ·  e^{-ik(R2-R1)}
        assert(ik == 0);
        for (auto cell : this->BvK_cells)
        {
            for (int it1 = 0;it1 < ucell.ntype;++it1)
                for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ++ia1)
                    for (int it2 = 0;it2 < ucell.ntype;++it2)
                        for (int ia2 = 0;ia2 < ucell.atoms[it2].na;++ia2)
                        {
                            int iat1 = ucell.itia2iat(it1, ia1);
                            int iat2 = ucell.itia2iat(it2, ia2);
                            auto& D2d = this->Ds_onebase[iat1][std::make_pair(iat2, cell)];
                            const int nw1 = aims_nbasis.empty() ? ucell.atoms[it1].nw : aims_nbasis[it1];
                            const int nw2 = aims_nbasis.empty() ? ucell.atoms[it2].nw : aims_nbasis[it2];
                            for (int iw1 = 0;iw1 < nw1;++iw1)
                                for (int iw2 = 0;iw2 < nw2;++iw2)
                                {
                                    const int iwt1 = ucell.itiaiw2iwt(it1, ia1, iw1);
                                    const int iwt2 = ucell.itiaiw2iwt(it2, ia2, iw2);
                                    if (this->pmat.in_this_processor(iwt1, iwt2))
                                        D2d(iw1, iw2) = this->psi_ks_full(ik, io, iwt1) * this->psi_ks_full(ik, nocc + iv, iwt2);
                                }
                        }
        }
    }

    template<>
    void OperatorLREXX<std::complex<double>>::cal_DM_onebase(const int io, const int iv, const int ik) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
        // NOTICE: DM_onebase will be passed into `cal_energy` interface and conjugated by "zdotc". 
        // So the formula should be the same as RHS. instead of LHS of the A-matrix, 
        // i.e. c1v · conj(c2o) ·  e^{-ik(R2-R1)}
        for (auto cell : this->BvK_cells)
        {
            std::complex<double> frac = RI::Global_Func::convert<std::complex<double>>(std::exp(
                -ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (this->kv.kvec_c.at(ik) * (RI_Util::array3_to_Vector3(cell) * GlobalC::ucell.latvec))));
            for (int it1 = 0;it1 < ucell.ntype;++it1)
                for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ++ia1)
                    for (int it2 = 0;it2 < ucell.ntype;++it2)
                        for (int ia2 = 0;ia2 < ucell.atoms[it2].na;++ia2)
                        {
                            int iat1 = ucell.itia2iat(it1, ia1);
                            int iat2 = ucell.itia2iat(it2, ia2);
                            auto& D2d = this->Ds_onebase[iat1][std::make_pair(iat2, cell)];
                            const int nw1 = aims_nbasis.empty() ? ucell.atoms[it1].nw : aims_nbasis[it1];
                            const int nw2 = aims_nbasis.empty() ? ucell.atoms[it2].nw : aims_nbasis[it2];
                            for (int iw1 = 0;iw1 < nw1;++iw1)
                                for (int iw2 = 0;iw2 < nw2;++iw2)
                                {
                                    const int iwt1 = ucell.itiaiw2iwt(it1, ia1, iw1);
                                    const int iwt2 = ucell.itiaiw2iwt(it2, ia2, iw2);
                                    if (this->pmat.in_this_processor(iwt1, iwt2))
                                        D2d(iw1, iw2) = frac * std::conj(this->psi_ks_full(ik, io, iwt2)) * this->psi_ks_full(ik, nocc + iv, iwt1);
                                }
                        }
        }
    }

    template<typename T>
    void OperatorLREXX<T>::act(const int nbands, const int nbasis, const int npol, const T* psi_in, T* hpsi, const int ngk_ik, const bool is_first_node)const
    {
        ModuleBase::TITLE("OperatorLREXX", "act");
        // convert parallel info to LibRI interfaces
        std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(this->pmat);

        // suppose Cs，Vs, have already been calculated in the ion-step of ground state
        // and DM_trans has been calculated in hPsi() outside.

        // 1. set_Ds (once)
        // convert to vector<T*> for the interface of RI_2D_Comm::split_m2D_ktoR (interface will be unified to ct::Tensor)
        std::vector<std::vector<T>> DMk_trans_vector = this->DM_trans->get_DMK_vector();
        // assert(DMk_trans_vector.size() == nk);
        std::vector<const std::vector<T>*> DMk_trans_pointer(nk);
        for (int ik = 0;ik < nk;++ik) { DMk_trans_pointer[ik] = &DMk_trans_vector[ik]; }
        // if multi-k, DM_trans(TR=double) -> Ds_trans(TR=T=complex<double>)
        std::vector<std::map<TA, std::map<TAC, RI::Tensor<T>>>> Ds_trans =
            aims_nbasis.empty() ?
            RI_2D_Comm::split_m2D_ktoR<T>(this->kv, DMk_trans_pointer, this->pmat, 1)
            : RI_Benchmark::split_Ds(DMk_trans_vector, aims_nbasis, ucell); //0.5 will be multiplied
        // LR_Util::print_CV(Ds_trans[0], "Ds_trans in OperatorLREXX", 1e-10);
        // 2. cal_Hs
        auto lri = this->exx_lri.lock();

        // LR_Util::print_CV(Ds_trans[is], "Ds_trans in OperatorLREXX", 1e-10);
        lri->exx_lri.set_Ds(std::move(Ds_trans[0]), lri->info.dm_threshold);
        lri->exx_lri.cal_Hs();
        lri->Hexxs[0] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
            lri->mpi_comm, std::move(lri->exx_lri.Hs), std::get<0>(judge[0]), std::get<1>(judge[0]));
        lri->post_process_Hexx(lri->Hexxs[0]);

        // 3. set [AX]_iak = DM_onbase * Hexxs for each occ-virt pair and each k-point
        // caution: parrallel

        for (int io = 0;io < this->nocc;++io)
        {
            for (int iv = 0;iv < this->nvirt;++iv)
            {
                for (int ik = 0;ik < nk;++ik)
                {
                    const int xstart_bk = ik * pX.get_local_size();
                    this->cal_DM_onebase(io, iv, ik);       //set Ds_onebase for all e-h pairs (not only on this processor)
                    // LR_Util::print_CV(Ds_onebase, "Ds_onebase of occ " + std::to_string(io) + ", virtual " + std::to_string(iv) + " in OperatorLREXX", 1e-10);
                    const T& ene = 2 * alpha * //minus for exchange(but here plus is right, why?), 2 for Hartree to Ry
                        lri->exx_lri.post_2D.cal_energy(this->Ds_onebase, lri->Hexxs[0]);
                    if (this->pX.in_this_processor(iv, io))
                    {
                        hpsi[xstart_bk + this->pX.global2local_col(io) * this->pX.get_row_size() + this->pX.global2local_row(iv)] += ene;
                    }
                }
            }
        }

    }
    template class OperatorLREXX<double>;
    template class OperatorLREXX<std::complex<double>>;
}
#endif
#ifdef __EXX
#include "operator_lr_exx.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{
    template<typename T>
    void OperatorLREXX<T>::allocate_Ds_onebase()
    {
        ModuleBase::TITLE("OperatorLREXX", "allocate_Ds_onebase");
        this->Ds_onebase.resize(this->nspin);
        for (int is = 0;is < this->nspin;++is) {
            for (int iat1 = 0;iat1 < ucell.nat;++iat1) {
                for (int iat2 = 0;iat2 < ucell.nat;++iat2) {
                    for (auto cell : this->BvK_cells) {
                        this->Ds_onebase[is][iat1][std::make_pair(iat2, cell)] =
                        RI::Tensor<T>({ static_cast<size_t>(ucell.atoms[ucell.iat2it[iat1]].nw),  static_cast<size_t>(ucell.atoms[ucell.iat2it[iat2]].nw) });
}
}
}
}
    }

    template<>
    void OperatorLREXX<double>::cal_DM_onebase(const int io, const int iv, const int ik, const int is) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
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
                            auto& D2d = this->Ds_onebase[is][iat1][std::make_pair(iat2, cell)];
                            for (int iw1 = 0;iw1 < ucell.atoms[it1].nw;++iw1)
                                for (int iw2 = 0;iw2 < ucell.atoms[it2].nw;++iw2)
                                    D2d(iw1, iw2) = this->psi_ks_full(ik, io, ucell.itiaiw2iwt(it1, ia1, iw1)) * this->psi_ks_full(ik, iv, ucell.itiaiw2iwt(it2, ia2, iw2));
                        }
        }
    }

    template<>
    void OperatorLREXX<std::complex<double>>::cal_DM_onebase(const int io, const int iv, const int ik, const int is) const
    {
        ModuleBase::TITLE("OperatorLREXX", "cal_DM_onebase");
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
                            auto& D2d = this->Ds_onebase[is][iat1][std::make_pair(iat2, cell)];
                            for (int iw1 = 0;iw1 < ucell.atoms[it1].nw;++iw1)
                                for (int iw2 = 0;iw2 < ucell.atoms[it2].nw;++iw2)
                                    D2d(iw1, iw2) = frac * std::conj(this->psi_ks_full(ik, io, ucell.itiaiw2iwt(it1, ia1, iw1))) * this->psi_ks_full(ik, iv, ucell.itiaiw2iwt(it2, ia2, iw2));
                        }
        }
    }

    template<typename T>
    void OperatorLREXX<T>::act(const psi::Psi<T>& psi_in, psi::Psi<T>& psi_out, const int nbands) const
    {
        ModuleBase::TITLE("OperatorLREXX", "act");

        assert(nbands <= psi_in.get_nbands());
        const int& nks = this->kv.get_nks();
        psi::Psi<T> psi_in_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_in, nks, this->pX->get_local_size());
        psi::Psi<T> psi_out_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_out, nks, this->pX->get_local_size());

        // convert parallel info to LibRI interfaces
        std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(*this->pmat);
        for (int ib = 0;ib < nbands;++ib)
        {
            psi_out_bfirst.fix_b(ib);
            // suppose Cs，Vs, have already been calculated in the ion-step of ground state, 
            // DM_trans(k) and DM_trans(R) has already been calculated from psi_in in OperatorLRHxc::act

            // 1. set_Ds (once)
            // convert to vector<T*> for the interface of RI_2D_Comm::split_m2D_ktoR (interface will be unified to ct::Tensor)
            std::cout << "ib=" << ib << std::endl;
            std::vector<std::vector<T>> DMk_trans_vector = this->DM_trans[ib]->get_DMK_vector();
            assert(DMk_trans_vector.size() == nks);
            std::vector<const std::vector<T>*> DMk_trans_pointer(nks);
            for (int is = 0;is < nks;++is) { DMk_trans_pointer[is] = &DMk_trans_vector[is];
}
            // if multi-k, DM_trans(TR=double) -> Ds_trans(TR=T=complex<double>)
            std::vector<std::map<TA, std::map<TAC, RI::Tensor<T>>>> Ds_trans =
                RI_2D_Comm::split_m2D_ktoR<T>(this->kv, DMk_trans_pointer, *this->pmat);

            // 2. cal_Hs
            for (int is = 0;is < nks;++is)
            {
                this->exx_lri->exx_lri.set_Ds(std::move(Ds_trans[is]), this->exx_lri->info.dm_threshold);
                this->exx_lri->exx_lri.cal_Hs();
                this->exx_lri->Hexxs[is] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
                    this->exx_lri->mpi_comm, std::move(this->exx_lri->exx_lri.Hs), std::get<0>(judge[is]), std::get<1>(judge[is]));
                this->exx_lri->post_process_Hexx(this->exx_lri->Hexxs[is]);
            }

            // 3. set [AX]_iak = DM_onbase * Hexxs for each occ-virt pair and each k-point
            // caution: parrallel

            for (int io = 0;io < this->pX->get_col_size();++io)   // nocc for serial
            {
                for (int iv = 0;iv < this->pX->get_row_size();++iv)   // nvirt for serial
                {
                    for (int ik = 0;ik < nks;++ik)
                    {
                        for (int is = 0;is < this->nspin;++is)
                        {
                            this->cal_DM_onebase(this->pX->local2global_col(io), this->pX->local2global_row(iv), ik, is);       //set Ds_onebase
                            psi_out_bfirst(ik, io * this->pX->get_row_size() + iv) -=
                                0.5 * //minus for exchange, 0.5 for spin
                                GlobalC::exx_info.info_global.hybrid_alpha *
                                this->exx_lri->exx_lri.post_2D.cal_energy(this->Ds_onebase[is], this->exx_lri->Hexxs[is]);
                        }
                    }
                }
            }
        }
    }
    template class OperatorLREXX<double>;
    template class OperatorLREXX<std::complex<double>>;
}
#endif
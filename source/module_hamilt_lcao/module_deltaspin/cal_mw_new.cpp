#include "spin_constrain.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_basis/module_ao/ORB_gen_tables.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW_new()
{
    std::cout << "hello world" << std::endl;
    const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
    const Parallel_Orbitals* pv = this->ParaV;
    int nat = this->get_nat();
    int npol = this->get_npol();
    /*
    for (int ipol = 0; ipol < this->npol_; ++ipol)
    {
        hamilt::HContainer<double>* dmR_current = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)
            ->get_DM()->get_DMR_pointer(ipol+1);
    }
    */
    for (const auto& sc_elem1 : this->get_atomCounts())
    {
        int it1 = sc_elem1.first;
        int nat_it1 = sc_elem1.second;
        for (int ia1 = 0; ia1 < nat_it1; ia1++)
        {
            int iat1 = this->get_iat(it1, ia1);
            auto tau1 = ucell->get_tau(iat1);
            AdjacentAtomInfo& adjs = this->adjs_all[iat1];
            std::vector<std::unordered_map<int, std::vector<double>>> nlm_tot;
            nlm_tot.resize(adjs.adj_num + 1);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int it2 = adjs.ntype[ad];
                const int ia2 = adjs.natom[ad];
                const int iat2 = this->get_iat(it2, ia2);
                const ModuleBase::Vector3<double>& tau2 = adjs.adjacent_tau[ad];
                const Atom* atom = &(ucell->atoms[it2]);
                auto all_indexes = pv->get_indexes_row(iat2);
                auto col_indexes = pv->get_indexes_col(iat2);
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
                for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
                {
                    const int iw1 = all_indexes[iw1l] / npol;
                    std::vector<std::vector<double>> nlm;
#ifdef USE_NEW_TWO_CENTER
                    int L = atom->iw2l[ iw1 ];
                    int N = atom->iw2n[ iw1 ];
                    int m = atom->iw2m[ iw1 ];
                    int M = (m % 2 == 0) ? -m/2 : (m+1)/2;
                    ModuleBase::Vector3<double> dtau = tau1 - tau2;
                    uot.two_center_bundle->overlap_orb->snap(it2, L, N, M, it1, dtau * ucell->lat0, 0, nlm);
#else
                    uot.snap_psibeta_half(orb,
                                        this->ucell->infoNL,
                                        nlm,
                                        tau2,
                                        it2,
                                        atom->iw2l[iw1],
                                        atom->iw2m[iw1],
                                        atom->iw2n[iw1],
                                        tau1,
                                        it1,
                                        0);
#endif
                    nlm_tot[ad].insert({all_indexes[iw1l], nlm[0]});
                }
            }
        }
    }
}
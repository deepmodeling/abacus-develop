#include "spin_constrain.h"
#include "module_basis/module_ao/ORB_read.h"

template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_adjs_all(const UnitCell& ucell, Grid_Driver* GridD)
{
    std::cout << "set_adjs_all" << std::endl;
    int nat = this->get_nat();
    this->adjs_all.clear();
    this->adjs_all.reserve(nat);
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
    for (const auto& sc_elem1 : this->get_atomCounts())
    {
        int it1 = sc_elem1.first;
        int nat_it1 = sc_elem1.second;
        for (int ia1 = 0; ia1 < nat_it1; ia1++)
        {
            int iat1 = this->get_iat(it1, ia1);
            auto tau1 = ucell.get_tau(iat1);
            AdjacentAtomInfo adjs;
            GridD->Find_atom(ucell, tau1, it1, ia1, &adjs);
            std::vector<bool> is_adj(adjs.adj_num + 1, false);
            for (int ad = 0; ad < adjs.adj_num + 1; ad++)
            {
                const int it2 = adjs.ntype[ad];
                const int ia2 = adjs.natom[ad];
                const int iat2 = this->get_iat(it2, ia2);
                const ModuleBase::Vector3<double>& tau2 = adjs.adjacent_tau[ad];
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                if (ucell.cal_dtau(iat1, iat2, R_index).norm()*ucell.lat0 <
                 orb.Phi[it2].getRcut() + ucell.infoNL.Beta[it1].get_rcut_max())
                {
                    is_adj[ad] = true;
                }
            }
            filter_adjs(is_adj, adjs);
            this->adjs_all.push_back(adjs);
        }
    }
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
template class SpinConstrain<double, psi::DEVICE_CPU>;
#include "../hcontainer.h"
void construct_SR(const Parallel_Orbitals* paraV, const UnitCell& ucell)
{
    this->SR = new hamilt::HContainer<double>(paraV);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int iat1=0;iat1<ucell.nat;iat++){
        auto tau1 = ucell.get_tau(iat1); 
        int T1 = ucell.get_element(iat1);
        int I1 = ucell.get_element_atom(iat1);
        GlobalC::GridD.Find_atom(ucell, tau1, T1, I1);
#ifdef _OPENMP
        hamilt::HContainer<double> SR_thread(paraV);
#endif
        for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
        {
            int iat2 = GlobalC::GridD.getAdjacent(ad);
            auto tau2 = ucell.get_tau(iat2); 
            int T2 = ucell.get_element(iat2);
            int I2 = ucell.get_element_atom(iat2);
            auto R_index = GlobalC::GridD.getBox(ad);
            double olm[3] = {0,0,0};
            hamilt::AtomPair<double> tmp(iat1, iat2, paraV);
            double* data_pointer = tmp.get_HR_values(R_index.x, R_index_y, R_index_z).get_pointer();
            for(iw1=0;iw2<ucell.get_nw(T1);++jj)
                for(iw2=0;iw2<ucell.get_nw(T2);++jj){
                    GlobalC::UOT.snap_psipsi( ... olm);
                    *data_pointer++ = olm[0];
                }
#ifdef _OPENMP
            SR_thread->insert_pair(tmp);
#else
            this->SR->insert_pair(tmp);
#endif
        }
#ifdef _OPENMP
#pragma omp critical(append_SR)
        this->SR->add_HR(SR_thread);
#endif
    }   
}
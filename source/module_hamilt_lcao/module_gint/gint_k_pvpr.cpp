#include "gint_k.h"
#include "grid_technique.h"
#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#ifdef __MPI
#include <mpi.h>
#endif

// transfer_pvpR, NSPIN = 1 or 2
void Gint_k::transfer_pvpR(hamilt::HContainer<double>* hR, const UnitCell* ucell, const Grid_Driver* gd)
{
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");

    for (int iap = 0; iap < this->hRGint->size_atom_pairs(); iap++)
    {
        auto& ap = this->hRGint->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = this->hRGint->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = this->hRGint->find_pair(iat1, iat2);
#ifdef __DEBUG
            assert(upper_ap != nullptr);
#endif
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {   
                auto R_index = ap.get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(-R_index);
                auto lower_mat = lower_ap->find_matrix(R_index);
                for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                    {
                        lower_mat->get_value(icol, irow) = upper_ap->get_value(irow, icol);
                    }
                }
            }
        }
    }
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*this->hRGint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*this->hRGint, hR);
    }
#else
    hR->add(*this->hRGint);
#endif
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    return;
}

// transfer_pvpR, NSPIN = 4
void Gint_k::transfer_pvpR(hamilt::HContainer<std::complex<double>>* hR,
                           const UnitCell* ucell_in,
                           const Grid_Driver* gd)
{
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");

    //this->hRGintCd->set_zero();
    
#ifdef __MPI
    int mg = hR->get_paraV()->get_global_row_size()/2;
    int ng = hR->get_paraV()->get_global_col_size()/2;
    int nb = hR->get_paraV()->get_block_size()/2;
    int blacs_ctxt = hR->get_paraV()->blacs_ctxt;
    int *iat2iwt = new int[ucell_in->nat];
    for (int iat = 0; iat < ucell_in->nat; iat++) {
        iat2iwt[iat] = ucell_in->get_iat2iwt()[iat]/2;
    }
    Parallel_Orbitals *pv = new Parallel_Orbitals();
    pv->set(mg, ng, nb, blacs_ctxt);
    pv->set_atomic_trace(iat2iwt, ucell_in->nat, mg);
    auto ijr_info = hR->get_ijr_info();
    

    hamilt::HContainer<double>* hR_tmp = new hamilt::HContainer<double>(pv, nullptr, &ijr_info);
    ModuleBase::Memory::record("Gint::hRGintCd", hR_tmp->get_memory_size());
    // std::cout<<"test1"<<std::endl;

    // std::cout<<hRGint_tmp[0]->size_atom_pairs()<<std::endl;
    // std::cout<<hR_tmp->size_atom_pairs()<<std::endl;
    // for(int i =0; i < hRGint_tmp[0]->size_atom_pairs(); i++){
    //     std::cout<<"hRGint_tmp"<<hRGint_tmp[0]->get_atom_pair(i).get_row_size()<<' '<<hRGint_tmp[0]->get_atom_pair(i).get_row_size()<<std::endl;
    //     std::cout<<"hR_tmp"<<hR_tmp->get_atom_pair(i).get_row_size()<<' '<<hR_tmp->get_atom_pair(i).get_row_size()<<std::endl;
    //     std::cout<<"pv1:" << pv->get_indexes_col(0).size()<< ' '<<pv->get_indexes_row(0).size()<<std::endl;
    //     std::cout<<"pv2:" << pv->get_indexes_col().size()<< ' '<<pv->get_indexes_row().size()<<std::endl;
    // }



    for (int is = 0; is < 4; is++){
        hR_tmp->set_zero();
        //std::cout<<"is: "<<is<<std::endl;
        hamilt::transferSerials2Parallels( *(this->hRGint_tmp[is]), hR_tmp);
        for (int iap = 0; iap < hR->size_atom_pairs(); iap++)
        {
            //std::cout<<"iap: "<<iap<<std::endl;
            auto* ap = &hR->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            const hamilt::AtomPair<double>* ap_nspin = nullptr;
            if (iat1 <= iat2)
            {
                hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
                hamilt::AtomPair<std::complex<double>>* lower_ap = hR->find_pair(iat2, iat1);
                switch (is)
                {
                case 0:
                    ap_nspin = hR_tmp->find_pair(iat1, iat2);
                    break;
                case 3:
                    ap_nspin = hR_tmp->find_pair(iat1, iat2);
                    break;
                }
                if(ap_nspin == nullptr) break;
                //const hamilt::AtomPair<double>* ap_nspin_0 = this->hR_tmp[0]->find_pair(iat1, iat2);
                //const hamilt::AtomPair<double>* ap_nspin_3 = this->hRGint_tmp[3]->find_pair(iat1, iat2);
                for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
                {   
                    const auto R_index = upper_ap->get_R_index(ir);
                    auto upper_mat = upper_ap->find_matrix(R_index);
                    //auto mat_nspin_0 = ap_nspin_0->find_matrix(R_index);
                    //auto mat_nspin_3 = ap_nspin_3->find_matrix(R_index);
                    auto mat_nspin = ap_nspin->find_matrix(R_index);

                    // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                    for (int irow = 0; irow < mat_nspin->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin->get_col_size(); ++icol)
                        {
                            //upper_mat->get_value(2*irow, 2*icol) = mat_nspin_0->get_value(irow, icol) + mat_nspin_3->get_value(irow, icol);
                            //upper_mat->get_value(2*irow+1, 2*icol+1) = mat_nspin_0->get_value(irow, icol) - mat_nspin_3->get_value(irow, icol);
                            switch (is)
                            {
                            case 0:
                                upper_mat->get_value(2*irow, 2*icol) = mat_nspin->get_value(irow, icol);
                                upper_mat->get_value(2*irow+1, 2*icol+1) = mat_nspin->get_value(irow, icol);
                                break;
                            case 3:
                                upper_mat->get_value(2*irow, 2*icol) += mat_nspin->get_value(irow, icol);
                                upper_mat->get_value(2*irow+1, 2*icol+1) -= mat_nspin->get_value(irow, icol);
                                break;
                            }
                        }
                    }

                    if (PARAM.globalv.domag)
                    {
                        const hamilt::AtomPair<double>* ap_nspin = nullptr;
                        switch (is)
                        {
                        case 1:
                            ap_nspin = hR_tmp->find_pair(iat1, iat2);
                            break;
                        case 2:
                            ap_nspin = hR_tmp->find_pair(iat1, iat2);
                            break;
                        }
                        //const hamilt::AtomPair<double>* ap_nspin_1 = this->hRGint_tmp[1]->find_pair(iat1, iat2);
                        //const hamilt::AtomPair<double>* ap_nspin_2 = this->hRGint_tmp[2]->find_pair(iat1, iat2);
                        //const auto mat_nspin_1 = ap_nspin_1->find_matrix(R_index);
                        //const auto mat_nspin_2 = ap_nspin_2->find_matrix(R_index);
                        const auto mat_nspin = ap_nspin->find_matrix(R_index);
                        for (int irow = 0; irow < mat_nspin->get_row_size(); ++irow)
                        {
                            for (int icol = 0; icol < mat_nspin->get_col_size(); ++icol)
                            {
                                switch(is)
                                {
                                    case 1:
                                        upper_mat->get_value(2*irow, 2*icol+1) = mat_nspin->get_value(irow, icol);
                                        upper_mat->get_value(2*irow+1, 2*icol) = mat_nspin->get_value(irow, icol);
                                        break;
                                    case 2:
                                        upper_mat->get_value(2*irow, 2*icol+1) += std::complex<double>(0.0, 1.0) * mat_nspin->get_value(irow, icol);
                                        upper_mat->get_value(2*irow+1, 2*icol) -= std::complex<double>(0.0, 1.0) * mat_nspin->get_value(irow, icol);
                                        break;
                                }
                                //upper_mat->get_value(2*irow, 2*icol+1) = mat_nspin->get_value(irow, icol) +  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                                //upper_mat->get_value(2*irow+1, 2*icol) = mat_nspin->get_value(irow, icol) -  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                            }
                        }
                    }
                    
                    // fill the lower triangle matrix
                    if(is == 3){
                        if (iat1 < iat2)
                        {
                            auto lower_mat = lower_ap->find_matrix(-R_index);
                            for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                            {
                                for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                                {
                                    lower_mat->get_value(icol, irow) = conj(upper_mat->get_value(irow, icol));
                                }
                            }
                        }
                    }
                        // std::cout<<"end"<<std::endl;
                }
            }
        }
        
    }
    delete[] iat2iwt;
    delete pv;
    delete hR_tmp;
#else

#endif

    // ===================================
    // transfer HR from Gint to Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>
    // ===================================
// #ifdef __MPI
//     int size;
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
//     if (size == 1)
//     {
//         hR->add(*this->hRGintCd);
//     }
//     else
//     {
//         hamilt::transferSerials2Parallels<std::complex<double>>(*this->hRGintCd, hR);
//     }
// #else
//     hR->add(*this->hRGintCd);
// #endif
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    return;
}

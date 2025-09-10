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

    int mg = hR->get_paraV()->get_global_row_size()/2;
    int ng = hR->get_paraV()->get_global_col_size()/2;
    int nb = hR->get_paraV()->get_block_size()/2;
#ifdef __MPI
    int blacs_ctxt = hR->get_paraV()->blacs_ctxt;
    int *iat2iwt = new int[ucell_in->nat];
    for (int iat = 0; iat < ucell_in->nat; iat++) {
        iat2iwt[iat] = ucell_in->get_iat2iwt()[iat]/2;
    }
    Parallel_Orbitals *pv = new Parallel_Orbitals();
    pv->set(mg, ng, nb, blacs_ctxt);
    pv->set_atomic_trace(iat2iwt, ucell_in->nat, mg);
    auto ijr_info = hR->get_ijr_info();

    this->hR_tmp = new hamilt::HContainer<std::complex<double>>(pv, nullptr, &ijr_info);
    ModuleBase::Memory::record("Gint::hRGintCd", this->hR_tmp->get_memory_size());

    //0,3;1,2;1,2;0,3
    std::vector<int> first = {0, 1, 1, 0};
    std::vector<int> second= {3, 2, 2, 3};
    for (int is = 0; is < 4; is++){
        if(!PARAM.globalv.domag && (is==1 || is==2)) continue;
        this->hR_tmp->set_zero();
        hamilt::HContainer<std::complex<double>>* hRGint_tmpCd = new hamilt::HContainer<std::complex<double>>(this->ucell->nat);
        hRGint_tmpCd->insert_ijrs(this->gridt->get_ijr_info(), *(this->ucell));
        hRGint_tmpCd->allocate(nullptr, true);
        hRGint_tmpCd->set_zero();
        for (int iap = 0; iap < hRGint_tmpCd->size_atom_pairs(); iap++)
        {
            //std::cout<<"iap: "<<iap<<std::endl;
            auto* ap = &hRGint_tmpCd->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            if (iat1 <= iat2)
            {
                hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
                hamilt::AtomPair<std::complex<double>>* lower_ap = hRGint_tmpCd->find_pair(iat2, iat1);
                const hamilt::AtomPair<double>* ap_nspin1 = this->hRGint_tmp[first[is]] ->find_pair(iat1, iat2);
                const hamilt::AtomPair<double>* ap_nspin2 = this->hRGint_tmp[second[is]] ->find_pair(iat1, iat2);
                for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
                {   
                    //std::cout<<"ir"<<ir<<std::endl;
                    const auto R_index = upper_ap->get_R_index(ir);
                    auto upper_mat = upper_ap->find_matrix(R_index);
                    auto mat_nspin1 = ap_nspin1->find_matrix(R_index);
                    auto mat_nspin2 = ap_nspin2->find_matrix(R_index);
                    // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                    for (int irow = 0; irow < mat_nspin1->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin1->get_col_size(); ++icol)
                        {
                            switch (is)
                            {
                            case 0:
                                upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) + mat_nspin2->get_value(irow, icol);
                                break;
                            case 1:
                                upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) + 
                                                                        std::complex<double>(0.0, 1.0) * mat_nspin2->get_value(irow, icol);
                                break;
                            case 2:
                                upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) - 
                                                                        std::complex<double>(0.0, 1.0) * mat_nspin2->get_value(irow, icol);
                                break;
                            case 3:
                                upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) - mat_nspin2->get_value(irow, icol);
                                break;
                            }
                        }
                    }
                    //fill the lower triangle matrix
                    if (is == 3 || is == 0){
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
                    if (is == 1 || is == 2){
                        if (iat1 < iat2)
                        {
                            auto lower_mat = lower_ap->find_matrix(-R_index);
                            for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                            {
                                for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                                {
                                    lower_mat->get_value(icol, irow) = upper_mat->get_value(irow, icol);
                                }
                            }
                        }
                    }
                }
            }
        }

        //std::cout<<"success"<<std::endl;

        hamilt::transferSerials2Parallels( *hRGint_tmpCd, this->hR_tmp);
        for (int iap = 0; iap < hR->size_atom_pairs(); iap++)
        {
            //std::cout<<"iap: "<<iap<<std::endl;
            auto* ap = &hR->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            auto* ap_nspin = this->hR_tmp ->find_pair(iat1, iat2);
            for (int ir = 0; ir < ap->get_R_size(); ir++)
            {   
                const auto R_index = ap->get_R_index(ir);
                auto upper_mat = ap->find_matrix(R_index);
                auto mat_nspin = ap_nspin->find_matrix(R_index);

                // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                for (int irow = 0; irow < mat_nspin->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < mat_nspin->get_col_size(); ++icol)
                    {
                        switch (is)
                        {
                        case 0:
                            upper_mat->get_value(2*irow, 2*icol) += mat_nspin->get_value(irow, icol);
                            break;
                        case 1:
                            upper_mat->get_value(2*irow, 2*icol+1) += mat_nspin->get_value(irow, icol);
                            break;
                        case 2:
                            upper_mat->get_value(2*irow+1, 2*icol) += mat_nspin->get_value(irow, icol);
                            break;
                        case 3:
                            upper_mat->get_value(2*irow+1, 2*icol+1) += mat_nspin->get_value(irow, icol);
                            break;
                        }
                    }
                }
            }
        }
    }
    delete[] iat2iwt;
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
    //  hR->add(*this->hRGintCd);
// #endif
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    return;
}

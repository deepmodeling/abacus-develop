#include "gint_common.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_parameter/parameter.h"

namespace ModuleGint
{

void compose_hRGint(std::shared_ptr<HContainer<double>> hRGint)
{
    for (int iap = 0; iap < hRGint->size_atom_pairs(); iap++)
    {
        auto& ap = hRGint->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = hRGint->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = hRGint->find_pair(iat1, iat2);
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
}

void compose_hRGint(std::vector<std::shared_ptr<HContainer<double>>> hRGint_part,
        std::shared_ptr<HContainer<std::complex<double>>> hRGint_full)
{
    for (int iap = 0; iap < hRGint_full->size_atom_pairs(); iap++)
    {
        auto* ap = &hRGint_full->get_atom_pair(iap);
        const int iat1 = ap->get_atom_i();
        const int iat2 = ap->get_atom_j();
        if (iat1 <= iat2)
        {
            hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
            hamilt::AtomPair<std::complex<double>>* lower_ap = hRGint_full->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* ap_nspin_0 = hRGint_part[0]->find_pair(iat1, iat2);
            const hamilt::AtomPair<double>* ap_nspin_3 = hRGint_part[3]->find_pair(iat1, iat2);
            for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
            {   
                const auto R_index = upper_ap->get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(R_index);
                auto mat_nspin_0 = ap_nspin_0->find_matrix(R_index);
                auto mat_nspin_3 = ap_nspin_3->find_matrix(R_index);

                // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                for (int irow = 0; irow < mat_nspin_0->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < mat_nspin_0->get_col_size(); ++icol)
                    {
                        upper_mat->get_value(2*irow, 2*icol) = mat_nspin_0->get_value(irow, icol) + mat_nspin_3->get_value(irow, icol);
                        upper_mat->get_value(2*irow+1, 2*icol+1) = mat_nspin_0->get_value(irow, icol) - mat_nspin_3->get_value(irow, icol);
                    }
                }

                if (PARAM.globalv.domag)
                {
                    const hamilt::AtomPair<double>* ap_nspin_1 = hRGint_part[1]->find_pair(iat1, iat2);
                    const hamilt::AtomPair<double>* ap_nspin_2 = hRGint_part[2]->find_pair(iat1, iat2);
                    const auto mat_nspin_1 = ap_nspin_1->find_matrix(R_index);
                    const auto mat_nspin_2 = ap_nspin_2->find_matrix(R_index);
                    for (int irow = 0; irow < mat_nspin_1->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin_1->get_col_size(); ++icol)
                        {
                            upper_mat->get_value(2*irow, 2*icol+1) = mat_nspin_1->get_value(irow, icol) +  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                            upper_mat->get_value(2*irow+1, 2*icol) = mat_nspin_1->get_value(irow, icol) -  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                        }
                    }
                }

                // fill the lower triangle matrix
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
        }
    }
}

template <typename T>
void transfer_hRGint_to_hR(std::shared_ptr<const HContainer<T>> hRGint, HContainer<T>* hR)
{
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*hRGint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*hRGint, hR);
    }
#else
    hR->add(*hRGint);
#endif
}

// gint_info should not have been a parameter, but it was added to initialize DMRGint_full
// In the future, we might try to remove the gint_info parameter
void transfer_DM_to_DMGint(
    std::shared_ptr<const GintInfo> gint_info,
    std::vector<HContainer<double>*> DM,
    std::vector<std::shared_ptr<HContainer<double>>> DMRGint)
{
    // To check whether input parameter DM2D has been initialized
#ifdef __DEBUG
    assert(PARAM.inp.nspin == DM.size()
           && "The size of DM should be equal to the number of spins!");
#endif

    if (PARAM.inp.nspin != 4)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
#ifdef __MPI
            hamilt::transferParallels2Serials(*DM[is], DMRGint[is].get());
#else
            DMRGint[is]->set_zero();
            DMRGint[is]->add(*DM[is]);
#endif
        }
    } else  // NSPIN=4 case
    {
#ifdef __MPI
        const int npol = 2;
        std::shared_ptr<HContainer<double>> DM_full = gint_info->get_hr<double>(npol);
        hamilt::transferParallels2Serials(*DM[0], DM_full.get());
#else
        HContainer<double>* DM_full = DM[0];
#endif
        std::vector<double*> tmp_pointer(4, nullptr);
        for (int iap = 0; iap < DM_full->size_atom_pairs(); iap++)
        {
            auto& ap = DM_full->get_atom_pair(iap);
            const int iat1 = ap.get_atom_i();
            const int iat2 = ap.get_atom_j();
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {
                const ModuleBase::Vector3<int> r_index = ap.get_R_index(ir);
                for (int is = 0; is < 4; is++)
                {
                    tmp_pointer[is] = 
                        DMRGint[is]->find_matrix(iat1, iat2, r_index)->get_pointer();
                }
                double* data_full = ap.get_pointer(ir);
                for (int irow = 0; irow < ap.get_row_size(); irow += 2)
                {
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2)
                    {
                        *(tmp_pointer[0])++ = data_full[icol];
                        *(tmp_pointer[1])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2)
                    {
                        *(tmp_pointer[2])++ = data_full[icol];
                        *(tmp_pointer[3])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                }
            }
        }
    }
}


template void transfer_hRGint_to_hR(std::shared_ptr<const HContainer<double>> hRGint, HContainer<double>* hR);
template void transfer_hRGint_to_hR(std::shared_ptr<const HContainer<std::complex<double>>> hRGint, HContainer<std::complex<double>>* hR);
}
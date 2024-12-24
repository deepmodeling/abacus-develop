#include "gint_common.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

namespace ModuleGint
{

void compose_hRGint(hamilt::HContainer<double>* hR)
{
    for (int iap = 0; iap < hR->size_atom_pairs(); iap++)
    {
        auto& ap = hR->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = hR->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = hR->find_pair(iat1, iat2);
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

template <typename T>
void transfer_hRGint_to_hR(const hamilt::HContainer<T>* hRGint, hamilt::HContainer<T>* hR)
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

template void transfer_hRGint_to_hR(const hamilt::HContainer<double>* hRGint, hamilt::HContainer<double>* hR);
template void transfer_hRGint_to_hR(const hamilt::HContainer<std::complex<double>>* hRGint, hamilt::HContainer<std::complex<double>>* hR);
}
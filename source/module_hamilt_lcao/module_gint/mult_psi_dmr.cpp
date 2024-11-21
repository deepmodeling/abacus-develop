#include "gint_tools.h"
#include "module_base/timer.h"
#include "module_base/ylm.h"

namespace Gint_Tools{

void mult_psi_DMR(
    const Grid_Technique& gt,
    const int bxyz,
    const int LD_pool,
    const int &grid_index,
    const int &na_grid,
    const int*const block_index,
    const int*const block_size,
    const bool*const*const cal_flag,
    const double*const*const psi,
    double*const*const psi_DMR,
    const hamilt::HContainer<double>*const DM,
    const bool if_symm)
{
    const UnitCell& ucell = *gt.ucell;

    // parameters for lapack subroutines
    constexpr char side = 'L';
    constexpr char uplo = 'U';
    const char trans = 'N';
    const double alpha = 1.0;
    const double beta = 1.0;
    const double alpha1 = if_symm ? 2.0 : 1.0;

    for (int ia1 = 0; ia1 < na_grid; ia1++)
    {
        const int bcell1 = gt.bcell_start[grid_index] + ia1;
        const int iat1 = gt.which_atom[bcell1];
        const int T1 = ucell.iat2it[iat1];
        const int I1 = ucell.iat2ia[iat1];

        //! get cell R1, this step is redundant in gamma_only case.
        const int id1 = gt.which_unitcell[bcell1];
        const int R1x = gt.ucell_index2x[id1];
        const int R1y = gt.ucell_index2y[id1];
        const int R1z = gt.ucell_index2z[id1];

        //! density
        if (if_symm)
        {
            //! ia2==ia1
            const auto tmp_matrix = DM->find_matrix(iat1, iat1, 0, 0, 0);
            
            //! maybe checking "tmp_matrix == nullptr" is not necessary
            if(tmp_matrix == nullptr)
            {
                continue;
            }
            
            const auto cal_info = Gint_Tools::cal_info(bxyz, ia1, ia1, cal_flag);
            const int ib_start = cal_info.first;
            const int ib_len = cal_info.second;
            
            if(ib_len == 0)
            {
                continue;
            }
            
            const auto tmp_matrix_ptr = tmp_matrix->get_pointer();
            const int idx1 = block_index[ia1];
            dsymm_(&side, &uplo, &block_size[ia1], &ib_len, &alpha, tmp_matrix_ptr, &block_size[ia1],
                    &psi[ib_start][idx1], &LD_pool, &beta, &psi_DMR[ib_start][idx1], &LD_pool);
        }

        //! get (j,beta,R2)
        const int start = if_symm ? ia1 + 1 : 0;

        for (int ia2 = start; ia2 < na_grid; ia2++)
        {
            const int bcell2 = gt.bcell_start[grid_index] + ia2;
            const int T2 = ucell.iat2it[gt.which_atom[bcell2]];
            const int iat2 = gt.which_atom[bcell2];
            const int id2 = gt.which_unitcell[bcell2];
 
            //! get cell R2, this step is redundant in gamma_only case.
            const int R2x = gt.ucell_index2x[id2];
            const int R2y = gt.ucell_index2y[id2];
            const int R2z = gt.ucell_index2z[id2];

            //! calculate the 'offset': R2 position relative
            //! to R1 atom, this step is redundant in gamma_only case.
            const int dRx = R1x - R2x;
            const int dRy = R1y - R2y;
            const int dRz = R1z - R2z;

            // get AtomPair
            const auto tmp_matrix = DM->find_matrix(iat1, iat2, dRx, dRy, dRz);
            if (tmp_matrix == nullptr)
            {
                continue;
            }
            const auto tmp_matrix_ptr = tmp_matrix->get_pointer();
            
            const auto cal_info = Gint_Tools::cal_info(bxyz, ia1, ia1, cal_flag);
            const int ib_start = cal_info.first;
            const int ib_len = cal_info.second;
            if(ib_len == 0)
            {
                continue;
            }
            const int idx1 = block_index[ia1];
            const int idx2 = block_index[ia2];
            
            dgemm_(&trans, &trans, &block_size[ia2], &ib_len, &block_size[ia1], &alpha1, tmp_matrix_ptr, &block_size[ia2],
                    &psi[ib_start][idx1], &LD_pool, &beta, &psi_DMR[ib_start][idx2], &LD_pool);

        }  // ia2
    } // ia1
}// End of mult_psi_DMR

}// End of Gint_Tools

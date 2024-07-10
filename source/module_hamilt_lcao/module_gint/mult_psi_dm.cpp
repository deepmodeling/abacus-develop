#include "gint_tools.h"

namespace Gint_Tools{
void mult_psi_DM(
    const Grid_Technique& gt, const int bxyz,
    const int na_grid, // how many atoms on this (i,j,k) grid
    const int LD_pool,
    const int* const block_iw,         // block_iw[na_grid],	index of wave functions for each block
    const int* const block_size,       // block_size[na_grid],	number of columns of a band
    const int* const block_index,      // block_index[na_grid+1], count total number of atomis orbitals
    const bool* const* const cal_flag, // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    const double* const* const psi,    // psir_vlbr3[bxyz][LD_pool]
    double** psi_DM, const double* const* const DM,
    const bool if_symm) // true: density, use dsymv; false: potential/transition density use dgemv
{
    constexpr char side = 'L', uplo = 'U';
    constexpr char transa = 'N', transb = 'N';
    constexpr double alpha_symm = 1, beta = 1;
    constexpr int inc = 1;
    double alpha_gemm = if_symm ? 2.0 : 1.0;

    for (int ia1 = 0; ia1 < na_grid; ia1++)
    {
        const int iw1_lo = block_iw[ia1];
        if (if_symm) // density
        {
            // ia1==ia2, diagonal part
            //  find the first ib and last ib for non-zeros cal_flag
            int first_ib = 0, last_ib = 0;
            for (int ib = 0; ib < bxyz; ++ib)
            {
                if (cal_flag[ib][ia1])
                {
                    first_ib = ib;
                    break;
                }
            }
            for (int ib = bxyz - 1; ib >= 0; --ib)
            {
                if (cal_flag[ib][ia1])
                {
                    last_ib = ib + 1;
                    break;
                }
            }
            const int ib_length = last_ib - first_ib;
            if (ib_length <= 0)
                continue;

            int cal_num = 0;
            for (int ib = first_ib; ib < last_ib; ++ib)
            {
                cal_num += cal_flag[ib][ia1];
            }
            // if enough cal_flag is nonzero
            if (cal_num > ib_length / 4)
            {
                dsymm_(&side, &uplo, &block_size[ia1], &ib_length, &alpha_symm, &DM[iw1_lo][iw1_lo], &gt.lgd,
                       &psi[first_ib][block_index[ia1]], &LD_pool, &beta, &psi_DM[first_ib][block_index[ia1]],
                       &LD_pool);
            }
            else
            {
                // int k=1;
                for (int ib = first_ib; ib < last_ib; ++ib)
                {
                    if (cal_flag[ib][ia1])
                    {
                        dsymv_(&uplo, &block_size[ia1], &alpha_symm, &DM[iw1_lo][iw1_lo], &gt.lgd,
                               &psi[ib][block_index[ia1]], &inc, &beta, &psi_DM[ib][block_index[ia1]], &inc);
                    }
                }
            }
        }

        int start = if_symm ? ia1 + 1 : 0;

        for (int ia2 = start; ia2 < na_grid; ia2++)
        {
            int first_ib = 0, last_ib = 0;
            for (int ib = 0; ib < bxyz; ++ib)
            {
                if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                {
                    first_ib = ib;
                    break;
                }
            }
            for (int ib = bxyz - 1; ib >= 0; --ib)
            {
                if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                {
                    last_ib = ib + 1;
                    break;
                }
            }
            const int ib_length = last_ib - first_ib;
            if (ib_length <= 0)
                continue;

            int cal_pair_num = 0;
            for (int ib = first_ib; ib < last_ib; ++ib)
            {
                cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
            }
            const int iw2_lo = block_iw[ia2];
            if (cal_pair_num > ib_length / 4)
            {
                dgemm_(&transa, &transb, &block_size[ia2], &ib_length, &block_size[ia1], &alpha_gemm,
                       &DM[iw1_lo][iw2_lo], &gt.lgd, &psi[first_ib][block_index[ia1]], &LD_pool, &beta,
                       &psi_DM[first_ib][block_index[ia2]], &LD_pool);
            }
            else
            {
                for (int ib = first_ib; ib < last_ib; ++ib)
                {
                    if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        dgemv_(&transa, &block_size[ia2], &block_size[ia1], &alpha_gemm, &DM[iw1_lo][iw2_lo], &gt.lgd,
                               &psi[ib][block_index[ia1]], &inc, &beta, &psi_DM[ib][block_index[ia2]], &inc);
                    }
                }
            }
        } // ia2
    }     // ia1
}
}
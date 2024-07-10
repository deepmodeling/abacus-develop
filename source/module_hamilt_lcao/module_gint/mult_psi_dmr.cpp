#include "gint_tools.h"

namespace Gint_Tools{
void mult_psi_DMR(const Grid_Technique& gt, const int bxyz, const int& grid_index, const int& na_grid,
                  const int* const block_index, const int* const block_size, bool** cal_flag, double** psi,
                  double** psi_DMR, const hamilt::HContainer<double>* DM, const bool if_symm)
{
    double *psi2, *psi2_dmr;
    int iwi, iww;
    const UnitCell& ucell = *gt.ucell;
    const int LD_pool = gt.max_atom * ucell.nwmax;

    // parameters for lapack subroutiens
    const char trans = 'N';
    const double alpha = 1.0, beta = 1.0;
    const int inc = 1;
    double alpha1;
    alpha1 = if_symm ? 2.0 : 1.0;

    for (int ia1 = 0; ia1 < na_grid; ia1++)
    {
        const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
        const int iat = gt.which_atom[mcell_index1];
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        Atom* atom1 = &ucell.atoms[T1];

        //~~~~~~~~~~~~~~~~
        // get cell R1.
        //~~~~~~~~~~~~~~~~
        const int id1 = gt.which_unitcell[mcell_index1];
        const int R1x = gt.ucell_index2x[id1];
        const int R1y = gt.ucell_index2y[id1];
        const int R1z = gt.ucell_index2z[id1];
        const double* tmp_matrix = DM->find_matrix(iat, iat, 0, 0, 0)->get_pointer();
        if (if_symm) // density
        {
            const int idx1 = block_index[ia1];
            std::vector<int> find_start = gt.find_R2[iat];
            std::vector<int>::const_iterator find_end = std::next(gt.find_R2[iat].begin(), gt.nad[iat]);
            // ia2==ia1
            int cal_num = 0;
            for (int ib = 0; ib < bxyz; ++ib)
            {
                if (cal_flag[ib][ia1])
                {
                    ++cal_num;
                }
            }

            int offset;
            if (cal_num > 0)
            {
                // find offset
                const int index = gt.cal_RindexAtom(0, 0, 0, iat);
                offset = -1;
                for (auto find = find_start.begin(); find < find_start.end(); find++)
                {
                    //--------------------------------------------------------------
                    // start positions of adjacent atom of 'iat'
                    //--------------------------------------------------------------
                    if (*find == index)
                    {
                        offset = find - find_start.begin(); // start positions of adjacent atom of 'iat'
                        break;
                    }
                }

                assert(offset != -1);
                assert(offset < gt.nad[iat]);
            }

            // const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];
            // for (int i = 0; i < block_size[ia1]; i++)
            //{
            //	for (int j = 0; j < block_size[ia1]; j++)
            //	{
            //		std::cout << i <<" "<<j<<" " << DMR[DM_start+i*block_size[ia1]+j]<<" " <<
            // tmp_matrix[i*block_size[ia1]+j] << std::endl;
            //	}
            // }
            // ModuleBase::WARNING_QUIT("cal_psi_dm","test");

            if (cal_num > bxyz / 4)
            {
                const int DM_start = gt.nlocstartg[iat] + gt.find_R2st[iat][offset];
                dgemm_(&trans, &trans, &block_size[ia1], &bxyz, &block_size[ia1], &alpha, tmp_matrix, &block_size[ia1],
                       &psi[0][idx1], &LD_pool, &beta, &psi_DMR[0][idx1], &LD_pool);
            }
            else if (cal_num > 0)
            {
                const int DM_start = gt.nlocstartg[iat] + gt.find_R2st[iat][offset];
                for (int ib = 0; ib < bxyz; ++ib)
                {
                    if (cal_flag[ib][ia1])
                    {
                        dgemv_(&trans, &block_size[ia1], &block_size[ia1], &alpha, tmp_matrix, &block_size[ia1],
                               &psi[ib][idx1], &inc, &beta, &psi_DMR[ib][idx1], &inc);
                    }
                }
            }
        }

        // get (j,beta,R2)
        int start = if_symm ? ia1 + 1 : 0;

        for (int ia2 = start; ia2 < na_grid; ia2++)
        {
            //---------------------------------------------
            // check if we need to calculate the big cell.
            //---------------------------------------------
            bool same_flag = false;
            for (int ib = 0; ib < gt.bxyz; ++ib)
            {
                if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                {
                    same_flag = true;
                    break;
                }
            }

            if (!same_flag)
                continue;

            const int bcell2 = gt.bcell_start[grid_index] + ia2;
            const int T2 = ucell.iat2it[gt.which_atom[bcell2]];
            const int iat2 = gt.which_atom[bcell2];

            Atom* atom2 = &ucell.atoms[T2];

            //---------------
            // get cell R2.
            //---------------
            const int id2 = gt.which_unitcell[bcell2];
            const int R2x = gt.ucell_index2x[id2];
            const int R2y = gt.ucell_index2y[id2];
            const int R2z = gt.ucell_index2z[id2];

            //------------------------------------------------
            // calculate the 'offset': R2 position relative
            // to R1 atom.
            //------------------------------------------------
            const int dRx = R1x - R2x;
            const int dRy = R1y - R2y;
            const int dRz = R1z - R2z;

            const int index = gt.cal_RindexAtom(dRx, dRy, dRz, iat2);
            // get AtomPair
            const double* tmp_matrix = DM->find_matrix(iat, iat2, dRx, dRy, dRz)->get_pointer();
            int offset = -1;

            std::vector<int> find_start(gt.find_R2[iat].begin(), gt.find_R2[iat].begin() + gt.nad[iat]);

            for (std::vector<int>::iterator find = find_start.begin(); find != find_start.end(); ++find)
            {
                if (*find == index)
                {
                    offset = std::distance(find_start.begin(), find);
                    break;
                }
            }
            if (offset == -1)
            {
                ModuleBase::WARNING_QUIT("gint_k", "mult_psi_DMR wrong");
            }
            assert(offset < gt.nad[iat]);

            //---------------------------------------------------------------
            // what I do above is to get 'offset' for atom std::pair (iat1, iat2)
            // if I want to simplify this searching for offset,
            // I should take advantage of gt.which_unitcell.
            //---------------------------------------------------------------

            int cal_num = 0;
            for (int ib = 0; ib < bxyz; ++ib)
            {
                if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    ++cal_num;
            }

            // const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];
            // for (int i = 0; i < block_size[ia1]; i++)
            //{
            //	for (int j = 0; j < block_size[ia1]; j++)
            //	{
            //		std::cout << i <<" "<<j<<" " << DMR[DM_start+i*block_size[ia1]+j]<<" " <<
            // tmp_matrix[i*block_size[ia1]+j] << std::endl;
            //	}
            // }
            // ModuleBase::WARNING_QUIT("cal_psi_dm","test");

            if (cal_num > bxyz / 4)
            {
                const int idx1 = block_index[ia1];
                const int idx2 = block_index[ia2];
                const int DM_start = gt.nlocstartg[iat] + gt.find_R2st[iat][offset];
                dgemm_(&trans, &trans, &block_size[ia2], &bxyz, &block_size[ia1], &alpha1, tmp_matrix, &block_size[ia2],
                       &psi[0][idx1], &LD_pool, &beta, &psi_DMR[0][idx2], &LD_pool);
            }
            else if (cal_num > 0)
            {
                const int idx1 = block_index[ia1];
                const int idx2 = block_index[ia2];
                const int DM_start = gt.nlocstartg[iat] + gt.find_R2st[iat][offset];

                for (int ib = 0; ib < bxyz; ++ib)
                {
                    if (cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        dgemv_(&trans, &block_size[ia2], &block_size[ia1], &alpha1, tmp_matrix, &block_size[ia2],
                               &psi[ib][idx1], &inc, &beta, &psi_DMR[ib][idx2], &inc);
                    }
                }
            } // cal_num
        }     // ia2
    }         // ia1
}
}
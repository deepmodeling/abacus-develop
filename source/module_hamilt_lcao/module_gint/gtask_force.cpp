#include <omp.h>

#include "gint_force.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
namespace GintKernel
{

/**
 * @brief Description of the function.
 *
 * Detailed description of the function.
 *
 * @param gridt The Grid_Technique object.
 * @param i The integer parameter i.
 * @param j The integer parameter j.
 * @param psiSizeMax The maximum size of psi.
 * @param max_size The maximum size.
 * @param nczp The nczp parameter.
 * @param vfactor The vfactor parameter.
 * @param vlocal_global_value The array of vlocal_global_value.
 * @param iat_per_nbz The array of iat_per_nbz.
 * @param input_dou The double array of psi_input.
 * @param psiInputInt The integer array of psi_input.
 * @param num_psir The array of num_psir.
 * @param lgd The lgd parameter.
 * @param psir_ylm_g The double array of psir_ylm_g.
 * @param psir_zeros_g The double array of psir_zeros_g.
 * @param dm_matrix_g The double array of dm_matrix_g.
 * @param mat_m The array of mat_m.
 * @param mat_n The array of mat_n.
 * @param mat_k The array of mat_k.
 * @param mat_lda The array of mat_lda.
 * @param mat_ldb The array of mat_ldb.
 * @param mat_ldc The array of mat_ldc.
 * @param mat_A The pointer to mat_A.
 * @param mat_B The pointer to mat_B.
 * @param mat_C The pointer to mat_C.
 * @param max_m The reference to max_m.
 * @param max_n The reference to max_n.
 * @param atom_pair_num The reference to atom_pair_num.
 */
void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int i,
                              const int j,
                              const int psiSizeMax,
                              const int max_size,
                              const int nczp,
                              const double vfactor,
                              double* rcut,
                              const double* vlocal_global_value,
                              int* iat_per_nbz,
                              const int lgd,
                              double* dm_matrix_g,
                              int& max_m,
                              int& max_n,
                              int& atom_pair_num,
                              SGridParameter& para)
{
    const int grid_index_ij = i * gridt.nby * gridt.nbzp + j * gridt.nbzp;
    const int nwmax = ucell.nwmax;
    bool* gpu_mat_cal_flag = new bool[max_size * gridt.nbzp];

    for (int i = 0; i < max_size * gridt.nbzp; i++)
    {
        gpu_mat_cal_flag[i] = false;
    }
    // psir generate
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int num_get_psi = 0;
        int grid_index = grid_index_ij + z_index;
        int num_psi_pos = psiSizeMax * z_index;
        int calc_flag_index = max_size * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int na_grid = gridt.how_many_atoms[grid_index];

        for (int id = 0; id < na_grid; id++)
        {
            int ib = 0;
            int mcell_index = bcell_start_index + id;
            int imcell = gridt.which_bigcell[mcell_index];
            int iat = gridt.which_atom[mcell_index];
            int it_temp = ucell.iat2it[iat];
            int start_ind_grid = gridt.start_ind[grid_index];

            for (int bx_index = 0; bx_index < gridt.bx; bx_index++)
            {
                for (int by_index = 0; by_index < gridt.by; by_index++)
                {
                    for (int bz_index = 0; bz_index < gridt.bz; bz_index++)
                    {
                        double dr_temp[3];
                        dr_temp[0] = gridt.meshcell_pos[ib][0]
                                     + gridt.meshball_positions[imcell][0]
                                     - gridt.tau_in_bigcell[iat][0];
                        dr_temp[1] = gridt.meshcell_pos[ib][1]
                                     + gridt.meshball_positions[imcell][1]
                                     - gridt.tau_in_bigcell[iat][1];
                        dr_temp[2] = gridt.meshcell_pos[ib][2]
                                     + gridt.meshball_positions[imcell][2]
                                     - gridt.tau_in_bigcell[iat][2];
                        /* compute distance in and allocate the paramter in
                         * z_index */
                        double distance = sqrt(dr_temp[0] * dr_temp[0]
                                               + dr_temp[1] * dr_temp[1]
                                               + dr_temp[2] * dr_temp[2]);
                        if (distance <= rcut[it_temp])
                        {
                            gpu_mat_cal_flag[calc_flag_index + id] = true;
                            int pos_temp_double = num_psi_pos + num_get_psi;
                            int pos_temp_int = pos_temp_double * 2;
                            pos_temp_double *= 5;
                            if (distance < 1.0E-9)
                                distance += 1.0E-9;
                            para.input_dou[pos_temp_double] = dr_temp[0];
                            para.input_dou[pos_temp_double + 1] = dr_temp[1];
                            para.input_dou[pos_temp_double + 2] = dr_temp[2];
                            para.input_dou[pos_temp_double + 3] = distance;
                            int vindex_global = bx_index * gridt.ncy * nczp
                                                + by_index * nczp + bz_index
                                                + start_ind_grid;
                            para.input_dou[pos_temp_double + 4]
                                = vlocal_global_value[vindex_global] * vfactor;

                            para.input_int[pos_temp_int] = it_temp;
                            para.input_int[pos_temp_int + 1]
                                = (z_index * gridt.bxyz + ib) * max_size * nwmax
                                  + id * nwmax;
                            iat_per_nbz[z_index * gridt.bxyz * max_size
                                        + ib * max_size + id]
                                = iat;
                            num_get_psi++;
                        }
                        ib++;
                    }
                }
            }
        }
        para.num_psir[z_index] = num_get_psi;
    }

    /* allocate the Multiplication of multinomial matrices */
    int tid = 0;
    max_m = 0;
    max_n = 0;

    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int calc_flag_index = max_size * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int bcell_start_psir = z_index * gridt.bxyz * max_size * nwmax;

        for (int atom1 = 0; atom1 < gridt.how_many_atoms[grid_index]; atom1++)
        {
            if (!gpu_mat_cal_flag[calc_flag_index + atom1])
            {
                continue;
            }
            int mcell_index1 = bcell_start_index + atom1;
            int iat1 = gridt.which_atom[mcell_index1];
            int it1 = ucell.iat2it[iat1];
            int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            int nw1 = ucell.atoms[it1].nw;

            for (int atom2 = 0; atom2 < gridt.how_many_atoms[grid_index];
                 atom2++)
            {
                if (!gpu_mat_cal_flag[calc_flag_index + atom2])
                {
                    continue;
                }
                int mcell_index2 = bcell_start_index + atom2;
                int iat2 = gridt.which_atom[mcell_index2];
                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                          ucell.iat2ia[iat2],
                                                          0)];
                int nw2 = ucell.atoms[it2].nw;

                int mat_A_idx = bcell_start_psir + atom2 * nwmax;
                int mat_B_idx = lgd * lo1 + lo2;
                int mat_C_idx = bcell_start_psir + atom1 * nwmax;
                para.atom_pair_A_m[tid] = gridt.bxyz;
                para.atom_pair_B_n[tid] = nw1;
                para.atom_pair_K[tid] = nw2;
                para.atom_pair_lda[tid] = nwmax * max_size;
                para.atom_pair_ldb[tid] = lgd;
                para.atom_pair_ldc[tid] = nwmax * max_size;
                para.matrix_A[tid] = para.psir_r_device + mat_A_idx;
                para.matrix_B[tid] = dm_matrix_g + mat_B_idx;
                para.matrix_C[tid] = para.psir_dm_device + mat_C_idx;

                if (para.atom_pair_A_m[tid] > max_m)
                {
                    max_m = para.atom_pair_A_m[tid];
                }

                if (para.atom_pair_B_n[tid] > max_n)
                {
                    max_n = para.atom_pair_B_n[tid];
                }

                tid++;
            }
        }
    }
    atom_pair_num = tid;

    delete[] gpu_mat_cal_flag;
}

void allocateDm(double* matrixHost,
                hamilt::HContainer<double>* dm,
                const Grid_Technique& gridt,
                const UnitCell& ucell)
{
    ModuleBase::GlobalFunc::ZEROS(matrixHost, gridt.lgd * gridt.lgd);
    for (int iatRow = 0; iatRow < ucell.nat; iatRow++)
    {
        for (int iatColumn = 0; iatColumn < ucell.nat; iatColumn++)
        {
            int indexTypeRow = ucell.iat2it[iatRow];
            int indexTypeColumn = ucell.iat2it[iatColumn];
            int localOrbitRow
                = gridt.trace_lo[ucell.itiaiw2iwt(indexTypeRow,
                                                  ucell.iat2ia[iatRow],
                                                  0)];
            int localOrbitColumn
                = gridt.trace_lo[ucell.itiaiw2iwt(indexTypeColumn,
                                                  ucell.iat2ia[iatColumn],
                                                  0)];
            hamilt::AtomPair<double>* tmpAtomPair
                = dm->find_pair(iatRow, iatColumn);
            int orbitIndex = 0;
            if (tmpAtomPair == NULL)
                continue;
            for (int orbitRow = 0; orbitRow < tmpAtomPair->get_row_size();
                 orbitRow++)
            {
                for (int orbitColumn = 0;
                     orbitColumn < tmpAtomPair->get_col_size();
                     orbitColumn++)
                {
                    matrixHost[(localOrbitRow + orbitRow) * gridt.lgd
                               + (localOrbitColumn + orbitColumn)]
                        = tmpAtomPair->get_pointer(0)[orbitIndex];
                    orbitIndex++;
                }
            }
        }
    }
    return;
}

} // namespace GintKernel

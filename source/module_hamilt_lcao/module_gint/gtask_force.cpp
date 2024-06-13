#include <omp.h>

#include "gint_force_gpu.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
namespace GintKernel
{

void gpu_task_generator_force(const Grid_Technique& gridt,
                              const UnitCell& ucell,
                              const int grid_index_ij,
                              const int max_atom_per_bcell,
                              const int max_atom,
                              const int nczp,
                              const double vfactor,
                              const double* rcut,
                              const double* vlocal_global_value,
                              double* psi_input_double,
                              int* psi_input_int,
                              int* atom_num_per_bcell,
                              int* iat_per_z,
                              int& atom_pair_num,
                              std::vector<bool>& gpu_mat_cal_flag)
{
    const int nwmax = ucell.nwmax;
    // psir generate
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int num_get_psi = 0;
        int grid_index = grid_index_ij + z_index;
        int num_psi_pos = max_atom_per_bcell * z_index;
        int calc_flag_index = max_atom * z_index;
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
                            {
                                distance += 1.0E-9;
                            }
                            psi_input_double[pos_temp_double] = dr_temp[0];
                            psi_input_double[pos_temp_double + 1] = dr_temp[1];
                            psi_input_double[pos_temp_double + 2] = dr_temp[2];
                            psi_input_double[pos_temp_double + 3] = distance;
                            int vindex_global = bx_index * gridt.ncy * nczp
                                                + by_index * nczp + bz_index
                                                + start_ind_grid;
                            psi_input_double[pos_temp_double + 4]
                                = vlocal_global_value[vindex_global] * vfactor;

                            psi_input_int[pos_temp_int] = it_temp;
                            psi_input_int[pos_temp_int + 1]
                                = (z_index * gridt.bxyz + ib) * max_atom * nwmax
                                  + id * nwmax;
                            iat_per_z[z_index * gridt.bxyz * max_atom
                                        + ib * max_atom + id]
                                = iat;
                            num_get_psi++;
                        }
                        ib++;
                    }
                }
            }
        }
        atom_num_per_bcell[z_index] = num_get_psi;
    }
}

/* allocate the Multiplication of multinomial matrices */
void alloc_mult_force(const Grid_Technique& gridt,
                      const UnitCell& ucell,
                      const int grid_index_ij,
                      const int max_atom,
                      double* const psi_g,
                      double* const psi_dm_g,
                      double* const dm_matrix_g,
                      int& max_m,
                      int& max_n,
                      int& atom_pair_num,
                      int* mat_m,
                      int* mat_n,
                      int* mat_k,
                      int* mat_lda,
                      int* mat_ldb,
                      int* mat_ldc,
                      double** mat_A,
                      double** mat_B,
                      double** mat_C,
                      std::vector<bool>& gpu_mat_cal_flag)
{
    int tid = 0;
    max_m = 0;
    max_n = 0;
    const int nwmax=ucell.nwmax;
    const int lgd = gridt.lgd;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int calc_flag_index = max_atom * z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int bcell_start_psir = z_index * gridt.bxyz * max_atom * nwmax;

        for (int atom1 = 0; atom1 < gridt.how_many_atoms[grid_index]; atom1++)
        {
            if (!gpu_mat_cal_flag[calc_flag_index + atom1])
            {
                continue;
            }
            const int mcell_index1 = bcell_start_index + atom1;
            int iat1 = gridt.which_atom[mcell_index1];
            int it1 = ucell.iat2it[iat1];
            int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            int nw1 = ucell.atoms[it1].nw;

            for (int atom2 = 0; atom2 < gridt.how_many_atoms[grid_index];atom2++)
            {
                if (!gpu_mat_cal_flag[calc_flag_index + atom2])
                {
                    continue;
                }
                const int mcell_index2 = bcell_start_index + atom2;
                int iat2 = gridt.which_atom[mcell_index2];
                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                            ucell.iat2ia[iat2],
                                                            0)];
                int nw2 = ucell.atoms[it2].nw;

                int mat_A_idx = bcell_start_psir + atom2 * nwmax;
                int mat_B_idx = lgd * lo1 + lo2;
                int mat_C_idx = bcell_start_psir + atom1 * nwmax;
                mat_m[tid] = gridt.bxyz;
                mat_n[tid] = nw1;
                mat_k[tid] = nw2;
                mat_lda[tid] = nwmax * max_atom;
                mat_ldb[tid] = lgd;
                mat_ldc[tid] = nwmax * max_atom;
                mat_A[tid] = psi_g + mat_A_idx;
                mat_B[tid] = dm_matrix_g + mat_B_idx;
                mat_C[tid] = psi_dm_g + mat_C_idx;

                if (mat_m[tid] > max_m)
                {
                    max_m = mat_m[tid];
                }

                if (mat_n[tid] > max_n)
                {
                    max_n = mat_n[tid];
                }

                tid++;
            }
        }
    }
    atom_pair_num = tid;
}
} // namespace GintKernel

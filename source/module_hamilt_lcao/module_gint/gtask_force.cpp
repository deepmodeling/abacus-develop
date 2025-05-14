#include <omp.h>

#include "gint_force_gpu.h"
#include "module_base/ylm.h"
#include "module_hamilt_lcao/module_gint/gint_tools.h"
#include "module_base/vector3.h"
namespace GintKernel
{

void gtask_force(const Grid_Technique& gridt,
                 const UnitCell& ucell,
                 const int grid_index_ij,
                 const int nczp,
                 const double vfactor,
                 const double* vlocal_global_value,
                 int& atoms_per_z,
                 int* atoms_num_info,
                 int* iat_on_nbz,
                 uint8_t* atoms_type,
                 double* dr_part,
                 double* vldr3)
{
    atoms_per_z = 0;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        const int grid_index = grid_index_ij + z_index;
        const int bcell_start_index = gridt.bcell_start[grid_index];
        const int na_grid = gridt.how_many_atoms[grid_index];
        atoms_num_info[z_index * 2] = na_grid;
        atoms_num_info[z_index * 2 + 1] = atoms_per_z;
        for (int id = 0; id < na_grid; id++)
        {
            const int mcell_index = bcell_start_index + id;
            const int imcell = gridt.which_bigcell[mcell_index];
            const int iat = gridt.which_atom[mcell_index];
            const int it_temp = ucell.iat2it[iat];

            dr_part[atoms_per_z * 3] = gridt.meshball_positions[imcell][0]
                                       - gridt.tau_in_bigcell[iat][0];
            dr_part[atoms_per_z * 3 + 1] = gridt.meshball_positions[imcell][1]
                                           - gridt.tau_in_bigcell[iat][1];
            dr_part[atoms_per_z * 3 + 2] = gridt.meshball_positions[imcell][2]
                                           - gridt.tau_in_bigcell[iat][2];
            atoms_type[atoms_per_z] = it_temp;
            iat_on_nbz[atoms_per_z] = iat;
            atoms_per_z++;
        }

        const int start_ind_grid = gridt.start_ind[grid_index];
        int id = z_index * gridt.bxyz;
        for (int bx_index = 0; bx_index < gridt.bx; bx_index++)
        {
            for (int by_index = 0; by_index < gridt.by; by_index++)
            {
                for (int bz_index = 0; bz_index < gridt.bz; bz_index++)
                {
                    int vindex_global = bx_index * gridt.ncy * nczp
                                        + by_index * nczp + bz_index
                                        + start_ind_grid;
                    vldr3[id]= vlocal_global_value[vindex_global] * vfactor;
                    id++;
                }
            }
        }
    }
}

void alloc_mult_force(const hamilt::HContainer<double>* dm,
                      const Grid_Technique& gridt,
                      const UnitCell& ucell,
                      const int grid_index_ij,
                      const int max_atom,
                      const int *atoms_num_info,
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
                      double** mat_C)
{
    int tid = 0;
    max_m = 0;
    max_n = 0;
    const int nwmax = ucell.nwmax;
    const int lgd = gridt.lgd;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        const int grid_index = grid_index_ij + z_index;
        const int bcell_start_index = gridt.bcell_start[grid_index];
        const int pre_atoms = atoms_num_info[z_index * 2 + 1];

        for (int atom1 = 0; atom1 < gridt.how_many_atoms[grid_index]; atom1++)
        {
            const int mcell_index1 = bcell_start_index + atom1;
            const int iat1 = gridt.which_atom[mcell_index1];
            const int uc1 = gridt.which_unitcell[mcell_index1];
            const ModuleBase::Vector3<int> r1 = gridt.get_ucell_coords(uc1); 
            const int it1 = ucell.iat2it[iat1];
            const int nw1 = ucell.atoms[it1].nw;

            for (int atom2 = 0; atom2 < gridt.how_many_atoms[grid_index];atom2++)
            {
                const int mcell_index2 = bcell_start_index + atom2;
                const int iat2 = gridt.which_atom[mcell_index2];
                const int uc2 = gridt.which_unitcell[mcell_index2];
                const ModuleBase::Vector3<int> r2 = gridt.get_ucell_coords(uc2);
                const int offset = dm->find_matrix_offset(iat1, iat2, r1-r2);
                if (offset == -1)
                {
                    continue;
                }
                const int it2 = ucell.iat2it[iat2];
                const int nw2 = ucell.atoms[it2].nw;

                const int mat_A_idx = (pre_atoms + atom2) * nwmax * gridt.bxyz;
                const int mat_C_idx = (pre_atoms + atom1) * nwmax * gridt.bxyz;
                mat_m[tid] = gridt.bxyz;
                mat_n[tid] = nw1;
                mat_k[tid] = nw2;
                mat_lda[tid] = nwmax;
                mat_ldb[tid] = nw2;
                mat_ldc[tid] = nwmax;
                mat_A[tid] = psi_g + mat_A_idx;
                mat_B[tid] = dm_matrix_g + offset;
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

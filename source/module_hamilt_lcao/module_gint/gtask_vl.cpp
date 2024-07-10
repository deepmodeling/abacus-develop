#include "gint_vl_gpu.h"

namespace GintKernel
{

void gtask_vlocal(const Grid_Technique& gridt,
                  const UnitCell& ucell,
                  const int grid_index_ij,
                  const int nczp,
                  const double vfactor,
                  const double* vlocal_global_value,
                  int& atoms_per_z,
                  int* atoms_num_info,
                  uint8_t* atoms_type,
                  double* dr_part,
                  double* vldr3)
{
    atoms_per_z = 0;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int bcell_start_index = gridt.bcell_start[grid_index];
        int na_grid = gridt.how_many_atoms[grid_index];
        atoms_num_info[2 * z_index] = na_grid;
        atoms_num_info[2 * z_index + 1] = atoms_per_z;
        for (int id = 0; id < na_grid; id++)
        {
            int mcell_index = bcell_start_index + id;
            int imcell = gridt.which_bigcell[mcell_index];
            int iat = gridt.which_atom[mcell_index];
            int it_temp = ucell.iat2it[iat];

            dr_part[atoms_per_z * 3] = gridt.meshball_positions[imcell][0]
                                       - gridt.tau_in_bigcell[iat][0];
            dr_part[atoms_per_z * 3 + 1] = gridt.meshball_positions[imcell][1]
                                           - gridt.tau_in_bigcell[iat][1];
            dr_part[atoms_per_z * 3 + 2] = gridt.meshball_positions[imcell][2]
                                           - gridt.tau_in_bigcell[iat][2];
            atoms_type[atoms_per_z] = it_temp;
            atoms_per_z++;
        }

        int start_ind_grid = gridt.start_ind[grid_index];
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

void alloc_mult_vlocal(const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        const int grid_index_ij,
                        const int max_atom,
                        double* psi,
                        double* psi_vldr3,
                        std::vector<Cuda_Mem_Wrapper<double>>& grid_vlocal_g,
                        int* mat_m,
                        int* mat_n,
                        int* mat_k,
                        int* mat_lda,
                        int* mat_ldb,
                        int* mat_ldc,
                        double** mat_A,
                        double** mat_B,
                        double** mat_C,
                        int& atom_pair_num,
                        int& max_m,
                        int& max_n)
{
    atom_pair_num = 0;
    max_m = 0;
    max_n = 0;
    const int nwmax = ucell.nwmax;
    for (int z_index = 0; z_index < gridt.nbzp; z_index++)
    {
        int grid_index = grid_index_ij + z_index;
        int atom_num = gridt.how_many_atoms[grid_index];
        int vldr3_index = z_index * max_atom * nwmax * gridt.bxyz;
        int bcell_start_index = gridt.bcell_start[grid_index];
        for (int atom1 = 0; atom1 < atom_num; atom1++)
        {
            int iat1 = gridt.which_atom[bcell_start_index + atom1];
            int it1 = ucell.iat2it[iat1];
            int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            for (int atom2 = 0; atom2 < atom_num; atom2++)
            {
                int iat2 = gridt.which_atom[bcell_start_index + atom2];
                int it2 = ucell.iat2it[iat2];
                int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                          ucell.iat2ia[iat2],
                                                          0)];
                if (lo1 <= lo2)
                {
                    int atom_pair_nw
                        = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                    if (grid_vlocal_g[iat1 * ucell.nat + iat2].get_device_pointer() == nullptr)
                    {
                        // Note that this situation occurs here because the
                        // logic in hcontainer and
                        //  grid integration is different.
                        //  In hcontainer, it is iat1<=iat2, and in grid
                        //  integral, it is lo1<=lo2. This is not entirely
                        //  equivalent in practice. We need to investigate
                        //  what's going on later.
                        //  TODO
                        continue;
                        // std::cout << "Error: GridVlocal did not malloc" <<
                        // std::endl;
                    }

                    int calc_index1 = vldr3_index + atom1 * nwmax * gridt.bxyz;
                    int calc_index2 = vldr3_index + atom2 * nwmax * gridt.bxyz;

                    mat_A[atom_pair_num]
                        = psi + calc_index1;
                    mat_B[atom_pair_num]
                        = psi_vldr3 + calc_index2;
                    mat_C[atom_pair_num]
                        = grid_vlocal_g[iat1 * ucell.nat + iat2].get_device_pointer();

                    mat_lda[atom_pair_num] = gridt.bxyz;
                    mat_ldb[atom_pair_num] = gridt.bxyz;
                    mat_ldc[atom_pair_num] = ucell.atoms[it2].nw;

                    mat_m[atom_pair_num] = ucell.atoms[it1].nw;
                    mat_n[atom_pair_num] = ucell.atoms[it2].nw;
                    mat_k[atom_pair_num] = gridt.bxyz;
                    
                    if (mat_m[atom_pair_num] > max_m)
                    {
                        max_m = mat_m[atom_pair_num];
                    }
                    if (mat_n[atom_pair_num] > max_n)
                    {
                        max_n = mat_n[atom_pair_num];
                    }
                    atom_pair_num++;
                }
            }
        }
    }
}

} // namespace GintKernel
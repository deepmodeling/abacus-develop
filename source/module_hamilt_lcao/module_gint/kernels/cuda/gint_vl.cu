#include "gint_vl.cuh"
#include "interp.cuh"
#include "cuda_tools.cuh"
#include "sph.cuh"
namespace GintKernel
{

__global__ void get_psi_and_vldr3(const double* const ylmcoef,
                                  const double delta_r,
                                  const int bxyz,
                                  const double nwmax,
                                  const double max_atom,
                                  const int* const ucell_atom_nwl,
                                  const bool* const atom_iw2_new,
                                  const int* const atom_iw2_ylm,
                                  const int* const atom_nw,
                                  const double* const rcut,
                                  const int nr_max,
                                  const double* const psi_u,
                                  const double* const mcell_pos,
                                  const double* const dr_part,
                                  const double* const vldr3,
                                  const uint8_t* const atoms_type,
                                  const int* const atoms_num_info,
                                  double* psi,
                                  double* psi_vldr3)
{
    const int bcell_id = blockIdx.x;
    const int num_atoms = atoms_num_info[2 * bcell_id];
    const int pre_atoms = atoms_num_info[2 * bcell_id + 1];
    const int mcell_id = blockIdx.y;
    const double vldr3_value = vldr3[bcell_id * bxyz + mcell_id];
    const double mcell_pos_x = mcell_pos[3 * mcell_id];
    const double mcell_pos_y = mcell_pos[3 * mcell_id + 1];
    const double mcell_pos_z = mcell_pos[3 * mcell_id + 2];

    for(int atom_id = threadIdx.x; atom_id < num_atoms; atom_id += blockDim.x)
    {
        const int dr_start = 3 * (pre_atoms + atom_id);
        const double dr_x = dr_part[dr_start] + mcell_pos_x;
        const double dr_y = dr_part[dr_start + 1] + mcell_pos_y;
        const double dr_z = dr_part[dr_start + 2] + mcell_pos_z;
        double dist = sqrt(dr_x * dr_x + dr_y * dr_y + dr_z * dr_z);
        const int atype = __ldg(atoms_type + pre_atoms + atom_id);
        if(dist < rcut[atype])
        {
            if (dist < 1.0E-9)
            {
                dist += 1.0E-9;
            }
            double dr[3] = {dr_x / dist, dr_y / dist, dr_z / dist};
            double ylma[49];
            const int nwl = __ldg(ucell_atom_nwl + atype);
            spherical_harmonics(dr, nwl, ylma, ylmcoef);
            int psi_idx = (bcell_id * max_atom + atom_id) * bxyz * nwmax + mcell_id;
            interp_vl(dist,
                      delta_r,
                      atype,
                      nwmax,
                      bxyz,
                      nr_max,
                      atom_nw,
                      atom_iw2_new,
                      psi_u,
                      ylma,
                      atom_iw2_ylm,
                      vldr3_value,
                      psi,
                      psi_vldr3,
                      psi_idx);
        }
    }
}

} // namespace GintKernel
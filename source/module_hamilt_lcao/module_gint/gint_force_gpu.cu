#include <omp.h>

#include "gint_force_gpu.h"
#include "kernels/cuda/cuda_tools.cuh"
#include "kernels/cuda/gint_force.cuh"
#include "module_base/ylm.h"
#include "gint_tools.h"

namespace GintKernel
{

// Function to calculate forces using GPU-accelerated gamma point Gint
/**
 * Function to calculate forces using GPU-accelerated gamma point Gint
 * @brief Calculate forces and stresses
 * @note The grid integration on the GPU is mainly divided into the following
 * steps:
 * 1. Use the CPU to divide the grid integration into subtasks.
 * 2. Copy the subtask information to the GPU.
 * 3. Calculate the matrix elements on the GPU.
 * 4. Perform matrix multiplication on the GPU.
 * 5. stress dot on the GPU.
 * 6. force dot on the GPU.
 * 7. Copy the results back to the host.
 */
void gint_fvl_gamma_gpu(hamilt::HContainer<double>* dm,
                          const double* vlocal,
                          double* force_in,
                          double* stress_in,
                          double dr,
                          double* rcut,
                          const int isforce,
                          const int isstress,
                          const Grid_Technique& gridt,
                          const UnitCell& ucell)
{   
    const int nbzp = gridt.nbzp;
    const int lgd = gridt.lgd;
    const int max_atom = gridt.max_atom;
    const int nwmax = ucell.nwmax;
    const int bxyz = gridt.bxyz;
    const int max_atom_per_bcell = max_atom * bxyz;
    const int max_atom_per_z = max_atom_per_bcell * nbzp;
    const int max_phi_per_z = max_atom_per_z * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;
    const double vfactor = ucell.omega / gridt.ncxyz;
    const int nczp = nbzp * gridt.bz;
    const int cuda_threads = 256;
    const int nat=ucell.nat;
    const int cuda_block
        = std::min(64, ceil_div(max_phi_per_z, cuda_threads));

    const int num_streams = gridt.nstreams;

    std::vector<cudaStream_t> streams(num_streams);
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }

    Cuda_Mem_Wrapper<double> psi_input_double(5 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> psi_input_int(2 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> atom_num_per_bcell(nbzp, num_streams, true);
    Cuda_Mem_Wrapper<int> start_idx_per_bcell(nbzp, num_streams, true);

    Cuda_Mem_Wrapper<double> psi(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> psi_dm(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> dpsi_dx(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> dpsi_dy(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> dpsi_dz(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dxx(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dxy(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dxz(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dyy(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dyz(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi_dzz(max_phi_per_z, num_streams, false);

    Cuda_Mem_Wrapper<double> gemm_alpha(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_m(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_n(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_k(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_lda(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_ldb(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_ldc(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_A(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_B(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_C(max_atompair_per_z, num_streams, true);

    Cuda_Mem_Wrapper<double> force(3 * nat, num_streams, true);
    Cuda_Mem_Wrapper<double> stress(6, num_streams, true);
    Cuda_Mem_Wrapper<int> iat_per_z(max_atom_per_z, num_streams, true);

    Cuda_Mem_Wrapper<double> dm_matrix(lgd * lgd, 1, true);
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < ucell.nat; iat2++)
        {
            const int it1 = ucell.iat2it[iat1];
            const int it2 = ucell.iat2it[iat2];
            const int lo1
                = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            const int lo2
                = gridt.trace_lo[ucell.itiaiw2iwt(it2, ucell.iat2ia[iat2], 0)];

            hamilt::AtomPair<double>* tmp_ap = dm->find_pair(iat1, iat2);
            int orb_index = 0;
            if (tmp_ap == NULL)
            {
                continue;
            }
            for (int orb_i = 0; orb_i < tmp_ap->get_row_size(); orb_i++)
            {
                for (int orb_j = 0; orb_j < tmp_ap->get_col_size(); orb_j++)
                {
                    dm_matrix.get_host_pointer()[(lo1 + orb_i) * lgd + (lo2 + orb_j)]
                        = tmp_ap->get_pointer(0)[orb_index];
                    orb_index++;
                }
            }
        }
    }
    dm_matrix.copy_host_to_device_sync();

    #pragma omp parallel for num_threads(num_streams) collapse(2)
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            const int sid = omp_get_thread_num();
            checkCuda(cudaStreamSynchronize(streams[sid]));

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;
            int atom_per_z = 0;
            const int grid_index_ij = i * gridt.nby * nbzp 
                                        + j * nbzp;

            std::vector<bool> gpu_mat_cal_flag(max_atom * nbzp, false);

            gtask_force(gridt,
                        ucell,
                        grid_index_ij,
                        max_atom_per_bcell,
                        max_atom,
                        nczp,
                        vfactor,
                        rcut,
                        vlocal,
                        psi_input_double.get_host_pointer(sid),
                        psi_input_int.get_host_pointer(sid),
                        atom_num_per_bcell.get_host_pointer(sid),
                        start_idx_per_bcell.get_host_pointer(sid),
                        iat_per_z.get_host_pointer(sid),
                        atom_per_z,
                        gpu_mat_cal_flag);
           
            alloc_mult_force(gridt,
                             ucell, 
                             grid_index_ij,
                             max_atom,
                             psi.get_device_pointer(sid),
                             psi_dm.get_device_pointer(sid),
                             dm_matrix.get_device_pointer(),
                             max_m,
                             max_n, 
                             atom_pair_num,
                             gemm_m.get_host_pointer(sid),
                             gemm_n.get_host_pointer(sid),
                             gemm_k.get_host_pointer(sid),
                             gemm_lda.get_host_pointer(sid),
                             gemm_ldb.get_host_pointer(sid),
                             gemm_ldc.get_host_pointer(sid),
                             gemm_A.get_host_pointer(sid),
                             gemm_B.get_host_pointer(sid),
                             gemm_C.get_host_pointer(sid),
                             gpu_mat_cal_flag);

            psi_input_double.copy_host_to_device_async(streams[sid], sid, 5 * atom_per_z);
            psi_input_int.copy_host_to_device_async(streams[sid], sid, 2 * atom_per_z);
            atom_num_per_bcell.copy_host_to_device_async(streams[sid], sid);
            start_idx_per_bcell.copy_host_to_device_async(streams[sid], sid);
            iat_per_z.copy_host_to_device_async(streams[sid], sid);
            gemm_m.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_n.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_k.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_lda.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_ldb.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_ldc.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_A.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_B.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_C.copy_host_to_device_async(streams[sid], sid, atom_pair_num);

            psi.memset_device_async(streams[sid], sid, 0);
            psi_dm.memset_device_async(streams[sid], sid, 0);
            dpsi_dx.memset_device_async(streams[sid], sid, 0);
            dpsi_dy.memset_device_async(streams[sid], sid, 0);
            dpsi_dz.memset_device_async(streams[sid], sid, 0);
            d2psi_dxx.memset_device_async(streams[sid], sid, 0);
            d2psi_dxy.memset_device_async(streams[sid], sid, 0);
            d2psi_dxz.memset_device_async(streams[sid], sid, 0);
            d2psi_dyy.memset_device_async(streams[sid], sid, 0);
            d2psi_dyz.memset_device_async(streams[sid], sid, 0);
            d2psi_dzz.memset_device_async(streams[sid], sid, 0);

            dim3 grid_psi(nbzp, 32);
            dim3 block_psi(64);
            get_psi_force<<<grid_psi,
                            block_psi,
                            0,
                            streams[sid]>>>(
                gridt.ylmcoef_g,
                dr,
                bxyz,
                nwmax,
                psi_input_double.get_device_pointer(sid),
                psi_input_int.get_device_pointer(sid),
                atom_num_per_bcell.get_device_pointer(sid),
                start_idx_per_bcell.get_device_pointer(sid),
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_l_g,
                gridt.atom_nw_g,
                gridt.nr_max,
                gridt.psi_u_g,
                psi.get_device_pointer(sid),
                dpsi_dx.get_device_pointer(sid),
                dpsi_dy.get_device_pointer(sid),
                dpsi_dz.get_device_pointer(sid),
                d2psi_dxx.get_device_pointer(sid),
                d2psi_dxy.get_device_pointer(sid),
                d2psi_dxz.get_device_pointer(sid),
                d2psi_dyy.get_device_pointer(sid),
                d2psi_dyz.get_device_pointer(sid),
                d2psi_dzz.get_device_pointer(sid));
            checkCudaLastError();

            gridt.fastest_matrix_mul(max_m,
                                     max_n,
                                     gemm_m.get_device_pointer(sid),
                                     gemm_n.get_device_pointer(sid),
                                     gemm_k.get_device_pointer(sid),
                                     gemm_A.get_device_pointer(sid),
                                     gemm_lda.get_device_pointer(sid),
                                     gemm_B.get_device_pointer(sid),
                                     gemm_ldb.get_device_pointer(sid),
                                     gemm_C.get_device_pointer(sid),
                                     gemm_ldc.get_device_pointer(sid),
                                     atom_pair_num,
                                     streams[sid],
                                     nullptr);

            if (isforce){
            const int block_size = std::min(256, ceil_div(nwmax, 32) * 32);
            dim3 grid_force(max_atom_per_z);
            dim3 block_force(block_size);
            dot_product_force<<<grid_force,
                                block_force,
                                block_size * 3 * sizeof(double),
                                streams[sid]>>>(
                                    dpsi_dx.get_device_pointer(sid),
                                    dpsi_dy.get_device_pointer(sid),
                                    dpsi_dz.get_device_pointer(sid),
                                    psi_dm.get_device_pointer(sid),
                                    force.get_device_pointer(sid),
                                    iat_per_z.get_device_pointer(sid),
                                    nwmax);
            checkCudaLastError();
            }

            if (isstress){ 
            dim3 grid_stress(cuda_block);
            dim3 block_stress(cuda_threads);
            dot_product_stress<<<grid_stress,
                                 block_stress,
                                 0,
                                 streams[sid]>>>(
                                d2psi_dxx.get_device_pointer(sid),
                                d2psi_dxy.get_device_pointer(sid),
                                d2psi_dxz.get_device_pointer(sid),
                                d2psi_dyy.get_device_pointer(sid),
                                d2psi_dyz.get_device_pointer(sid),
                                d2psi_dzz.get_device_pointer(sid),
                                psi_dm.get_device_pointer(sid),
                                stress.get_device_pointer(sid),
                                max_phi_per_z);
            checkCudaLastError();
            }
        }
    }

    for(int i = 0; i < num_streams; i++)
    {
        stress.copy_device_to_host_async(streams[i], i);
        force.copy_device_to_host_async(streams[i], i);
    }

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamSynchronize(streams[i]));
    }

    if (isstress){
        for (int i = 0; i < num_streams; i++)
        {
            const int offset = 6 * i;
            for (int j = 0; j < 6; j++)
            {
                stress_in[j] += stress.get_host_pointer()[offset + j];
            }
        }
    }
    if (isforce){
        for (int i = 0; i < num_streams; i++)
        {
            const int offset = 3 * i * nat;
            for (int j = 0; j < 3 * nat; j++)
            {
                force_in[j] += force.get_host_pointer()[offset + j];
            }
        }
    }

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamDestroy(streams[i]));
    }
}

} // namespace GintKernel

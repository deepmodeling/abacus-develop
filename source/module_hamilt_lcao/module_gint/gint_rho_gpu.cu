#include "kernels/cuda/cuda_tools.cuh"
#include "module_base/ylm.h"
#include "gint_rho_gpu.h"
#include "gint_tools.h"
#include "kernels/cuda/gint_rho.cuh"

#include <omp.h>

namespace GintKernel
{

void gint_gamma_rho_gpu(const hamilt::HContainer<double>* dm,
                        const double* ylmcoef_now,
                        const double dr,
                        const double* rcut,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        double* rho)
{
    const int nbzp = gridt.nbzp;
    const int nczp =nbzp * gridt.bz;
    const int num_mcell_on_proc = nczp * gridt.ncx * gridt.ncy;
    const int lgd = gridt.lgd;
    const int max_atom = gridt.max_atom;
    const int num_streams = gridt.nstreams;
    const int max_atom_per_bcell = max_atom * gridt.bxyz;
    const int max_atom_per_z = max_atom_per_bcell * nbzp;
    const int max_phi_per_z = max_atom_per_z * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;

    std::vector<cudaStream_t> streams(num_streams);
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }
    
    Cuda_Mem_Wrapper<double> psi_input_double(4 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> psi_input_int(2 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> atom_num_per_bcell(nbzp, num_streams, true);
    Cuda_Mem_Wrapper<int> start_idx_per_bcell(nbzp, num_streams, true);

    Cuda_Mem_Wrapper<double> psi(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> psi_dm(max_phi_per_z, num_streams, false);

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
    
    Cuda_Mem_Wrapper<double> rho_g(num_mcell_on_proc, 1, false);
    Cuda_Mem_Wrapper<double*> dot_product(nbzp * gridt.bxyz, num_streams, true);

    Cuda_Mem_Wrapper<double> dm_matrix(lgd * lgd, 1, true);
    // retrieve the density matrix on the host
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < ucell.nat; iat2++)
        {
            const int it1 = ucell.iat2it[iat1];
            const int it2 = ucell.iat2it[iat2];
            const int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(it1, ucell.iat2ia[iat1], 0)];
            const int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2, ucell.iat2ia[iat2], 0)];

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

// calculate the rho for every nbzp bigcells
#pragma omp parallel for num_threads(num_streams) collapse(2)
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            // get stream id
            const int sid = omp_get_thread_num();
            checkCuda(cudaStreamSynchronize(streams[sid]));

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;
            int atom_per_z = 0;
            const int grid_index_ij = i * gridt.nby * nbzp + j * nbzp;
            std::vector<bool> gpu_matrix_cal_flag(max_atom * nbzp, false);

            // generate GPU tasks, including the calculation of psir, matrix
            // multiplication, and dot product
            gtask_rho(gridt,
                      grid_index_ij,
                      gpu_matrix_cal_flag,
                      max_atom,
                      ucell,
                      rcut,
                      psi_input_double.get_host_pointer(sid),
                      psi_input_int.get_host_pointer(sid),
                      atom_num_per_bcell.get_host_pointer(sid),
                      start_idx_per_bcell.get_host_pointer(sid),
                      atom_per_z);
            
            alloc_mult_dot_rho(gridt,
                            ucell,
                            gpu_matrix_cal_flag,
                            grid_index_ij,
                            max_atom,
                            lgd,
                            nczp,
                            psi.get_device_pointer(sid),
                            psi_dm.get_device_pointer(sid),
                            dm_matrix.get_device_pointer(),
                            gemm_alpha.get_host_pointer(sid),
                            gemm_m.get_host_pointer(sid),
                            gemm_n.get_host_pointer(sid),
                            gemm_k.get_host_pointer(sid),
                            gemm_lda.get_host_pointer(sid),
                            gemm_ldb.get_host_pointer(sid),
                            gemm_ldc.get_host_pointer(sid),
                            gemm_A.get_host_pointer(sid),
                            gemm_B.get_host_pointer(sid),
                            gemm_C.get_host_pointer(sid),
                            max_m,
                            max_n,
                            atom_pair_num,
                            rho_g.get_device_pointer(),
                            dot_product.get_host_pointer(sid));
           
            psi_input_double.copy_host_to_device_async(streams[sid], sid, 4 * atom_per_z);
            psi_input_int.copy_host_to_device_async(streams[sid], sid, 2 * atom_per_z);
            atom_num_per_bcell.copy_host_to_device_async(streams[sid], sid);
            start_idx_per_bcell.copy_host_to_device_async(streams[sid], sid);
            gemm_alpha.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_m.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_n.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_k.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_lda.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_ldb.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_ldc.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_A.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_B.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            gemm_C.copy_host_to_device_async(streams[sid], sid, atom_pair_num);
            dot_product.copy_host_to_device_async(streams[sid], sid);
            
            psi.memset_device_async(streams[sid], sid, 0);
            psi_dm.memset_device_async(streams[sid], sid, 0);

            // Launching kernel to calculate psi
            dim3 grid_psi(nbzp, 8);
            dim3 block_psi(64);
            get_psi<<<grid_psi, block_psi, 0, streams[sid]>>>(
                gridt.ylmcoef_g,
                dr,
                gridt.bxyz,
                ucell.nwmax,
                psi_input_double.get_device_pointer(sid),
                psi_input_int.get_device_pointer(sid),
                atom_num_per_bcell.get_device_pointer(sid),
                start_idx_per_bcell.get_device_pointer(sid),
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_nw_g,
                gridt.nr_max,
                gridt.psi_u_g,
                psi.get_device_pointer(sid));
            checkCudaLastError();

            // Performing matrix multiplication alpha * mat_dm * mat_psir
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
                                     gemm_alpha.get_device_pointer(sid));
            checkCudaLastError();

            // Launching kernel to calculate dot product psir * psir_dm
            const int block_size = 128;
            dim3 block_dot(block_size);
            dim3 grid_dot(nbzp, gridt.bxyz);
            psir_dot<<<grid_dot, block_dot, sizeof(double) * block_size, streams[sid]>>>(
                gridt.bxyz,
                max_atom * ucell.nwmax,
                psi.get_device_pointer(sid),
                psi_dm.get_device_pointer(sid),
                dot_product.get_device_pointer(sid));
            checkCudaLastError();
        }
    }

    // Synchronizing streams
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamSynchronize(streams[i]));
    }

    // Copy rho from device to host
    checkCuda(cudaMemcpy(rho,
                         rho_g.get_device_pointer(),
                         num_mcell_on_proc * sizeof(double),
                         cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamDestroy(streams[i]));
    }
}
} // namespace GintKernel

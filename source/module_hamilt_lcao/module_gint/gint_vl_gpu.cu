#include <omp.h>

#include "kernels/cuda/cuda_tools.cuh"
#include "module_base/ylm.h"
#include "gint_tools.h"
#include "gint_vl_gpu.h"
#include "kernels/cuda/gint_vl.cuh"

namespace GintKernel
{

/**
 * Computes the gamma component of the VL (Vlocal) integral on the GPU.
 *
 * @param hRGint Pointer to the HContainer<double> object to store the computed
 * integrals.
 * @param lgd Dimension information for the computation results.
 * @param max_atom The maximum number of neighboring atoms for a grid point.
 * @param vfactor Related to volume. The scaling factor for the Vlocal
 * integrals.
 * @param vlocal Pointer to the Vlocal array.
 * @param ylmcoef_now Pointer to the Ylm coefficients array.
 * @param nczp The number of grid layers in the C direction.
 * @param nbxx The total number of grid points.
 * @param dr The grid spacing.
 * @param rcut Pointer to the cutoff radius array.
 * @param gridt The Grid_Technique object containing grid information.
 * @param ucell The UnitCell object containing unit cell information.
 *
 * @note The grid integration on the GPU is mainly divided into the following
 * steps:
 * 1. Use the CPU to divide the grid integration into subtasks.
 * 2. Copy the subtask information to the GPU.
 * 3. Calculate the matrix elements on the GPU.
 * 4. Perform matrix multiplication on the GPU.
 * 5. Copy the results back to the host.
 */
void gint_gamma_vl_gpu(hamilt::HContainer<double>* hRGint,
                       const double* vlocal,
                       const double* ylmcoef_now,
                       const double dr,
                       const double* rcut,
                       const Grid_Technique& gridt,
                       const UnitCell& ucell)
{
    const int nbzp = gridt.nbzp;
    const int num_streams = gridt.nstreams;
    const int lgd = gridt.lgd;
    const int max_atom = gridt.max_atom;
    const int max_atom_per_bcell = max_atom * gridt.bxyz;
    const int max_atom_per_z = max_atom_per_bcell * nbzp;
    const int max_phi_per_z = max_atom_per_z * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;
    const double vfactor = ucell.omega / gridt.ncxyz;
    const int nczp = nbzp * gridt.bz;
    const int nbxx = gridt.nbxx;  // total number of grid points
    std::vector<cudaStream_t> streams(num_streams);

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }

    std::vector<Cuda_Mem_Wrapper<double>> grid_vlocal_g(ucell.nat * ucell.nat);
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        for (int iat2 = 0; iat2 < ucell.nat; iat2++)
        {
            const int it1 = ucell.iat2it[iat1];
            const int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(it1,
                                                        ucell.iat2ia[iat1],
                                                        0)];

            const int it2 = ucell.iat2it[iat2];
            const int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                        ucell.iat2ia[iat2],
                                                        0)];

            if (lo1 <= lo2)
            {
                const hamilt::AtomPair<double>* tmp_ap
                    = hRGint->find_pair(iat1, iat2);
                if (tmp_ap == nullptr)
                {
                    continue;
                }
                const int atom_pair_nw
                    = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                grid_vlocal_g[iat1 * ucell.nat + iat2] = 
                    Cuda_Mem_Wrapper<double>(atom_pair_nw, 1, false);
                grid_vlocal_g[iat1 * ucell.nat + iat2].memset_device_sync();
            }
        }
    }

    Cuda_Mem_Wrapper<double> psi_input_double(5 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> psi_input_int(2 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> atom_num_per_bcell(nbzp, num_streams, true);

    Cuda_Mem_Wrapper<double> psi(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> psi_vldr3(max_phi_per_z, num_streams, false);

    Cuda_Mem_Wrapper<int> gemm_m(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_n(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_k(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_lda(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_ldb(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> gemm_ldc(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_A(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_B(max_atompair_per_z, num_streams, true);
    Cuda_Mem_Wrapper<double*> gemm_C(max_atompair_per_z, num_streams, true);
    
#pragma omp parallel for num_threads(num_streams) collapse(2)
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            const int sid = omp_get_thread_num();
            checkCuda(cudaStreamSynchronize(streams[sid]));

            int atom_pair_num = 0;
            int max_m = 0;
            int max_n = 0;
            const int grid_index_ij = i * gridt.nby * gridt.nbzp + j * gridt.nbzp;
            std::vector<bool> gpu_matrix_calc_flag(max_atom * nbzp, false);
            gtask_vlocal(gridt,
                         rcut,
                         ucell,
                         gpu_matrix_calc_flag,
                         grid_index_ij,
                         max_atom,
                         nczp,
                         vfactor,
                         vlocal,
                         psi_input_double.get_host_pointer(sid),
                         psi_input_int.get_host_pointer(sid),
                         atom_num_per_bcell.get_host_pointer(sid));
        
            alloc_mult_vlocal(gridt,
                                ucell,
                                gpu_matrix_calc_flag,
                                grid_index_ij,
                                max_atom,
                                psi.get_device_pointer(sid),
                                psi_vldr3.get_device_pointer(sid),
                                grid_vlocal_g,
                                gemm_m.get_host_pointer(sid),
                                gemm_n.get_host_pointer(sid),
                                gemm_k.get_host_pointer(sid),
                                gemm_lda.get_host_pointer(sid),
                                gemm_ldb.get_host_pointer(sid),
                                gemm_ldc.get_host_pointer(sid),
                                gemm_A.get_host_pointer(sid),
                                gemm_B.get_host_pointer(sid),
                                gemm_C.get_host_pointer(sid),
                                atom_pair_num,
                                max_m,
                                max_n);

            psi_input_double.copy_host_to_device_async(streams[sid], sid);
            psi_input_int.copy_host_to_device_async(streams[sid], sid);
            atom_num_per_bcell.copy_host_to_device_async(streams[sid], sid);
            gemm_m.copy_host_to_device_async(streams[sid], sid);
            gemm_n.copy_host_to_device_async(streams[sid], sid);
            gemm_k.copy_host_to_device_async(streams[sid], sid);
            gemm_lda.copy_host_to_device_async(streams[sid], sid);
            gemm_ldb.copy_host_to_device_async(streams[sid], sid);
            gemm_ldc.copy_host_to_device_async(streams[sid], sid);
            gemm_A.copy_host_to_device_async(streams[sid], sid);
            gemm_B.copy_host_to_device_async(streams[sid], sid);
            gemm_C.copy_host_to_device_async(streams[sid], sid);

            psi.memset_device_async(streams[sid], sid, 0);
            psi_vldr3.memset_device_async(streams[sid], sid, 0);

            dim3 grid_psi(nbzp, 8);
            dim3 block_psi(64);
            get_psi_and_vldr3<<<grid_psi,
                                block_psi,
                                0,
                                streams[sid]>>>(
                gridt.ylmcoef_g,
                dr,
                gridt.bxyz,
                ucell.nwmax,
                psi_input_double.get_device_pointer(sid),
                psi_input_int.get_device_pointer(sid),
                atom_num_per_bcell.get_device_pointer(sid),
                max_atom_per_bcell,
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_nw_g,
                gridt.nr_max,
                gridt.psi_u_g,
                psi.get_device_pointer(sid),
                psi_vldr3.get_device_pointer(sid));
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
            checkCudaLastError();
        }
    }

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamSynchronize(streams[i]));
    }
    {
        int iter_num = 0;
        for (int iat1 = 0; iat1 < ucell.nat; iat1++)
        {
            for (int iat2 = 0; iat2 < ucell.nat; iat2++)
            {
                const int sid = iter_num % num_streams;
                const int it1 = ucell.iat2it[iat1];
                const int lo1 = gridt.trace_lo[ucell.itiaiw2iwt(it1,
                                                          ucell.iat2ia[iat1],
                                                          0)];

                const int it2 = ucell.iat2it[iat2];
                const int lo2 = gridt.trace_lo[ucell.itiaiw2iwt(it2,
                                                          ucell.iat2ia[iat2],
                                                          0)];
                if (lo1 <= lo2)
                {
                    const int atom_pair_nw
                        = ucell.atoms[it1].nw * ucell.atoms[it2].nw;
                    hamilt::AtomPair<double>* tmp_ap
                        = hRGint->find_pair(iat1, iat2);
                    if (tmp_ap == nullptr)
                    {
                        continue;
                    }
                    checkCuda(cudaMemcpyAsync(
                        tmp_ap->get_pointer(0),
                        grid_vlocal_g[iat1 * ucell.nat + iat2].get_device_pointer(),
                        atom_pair_nw * sizeof(double),
                        cudaMemcpyDeviceToHost,
                        streams[sid]));
                    iter_num++;
                }
            }
        }
    }
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamSynchronize(streams[i]));
    }
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamDestroy(streams[i]));
    }
}

} // namespace GintKernel
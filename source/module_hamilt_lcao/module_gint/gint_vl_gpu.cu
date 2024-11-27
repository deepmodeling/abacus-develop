#ifdef _OPENMP
#include <omp.h>
#endif

#include "kernels/cuda/cuda_tools.cuh"
#include "module_base/ylm.h"
#include "gint_vl_gpu.h"
#include "kernels/cuda/gint_vl.cuh"

namespace GintKernel
{

/**
 * Computes the gamma component of the VL (Vlocal) integral on the GPU.
 *
 * @note The grid integration on the GPU is mainly divided into the following
 * steps:
 * 1. Use the CPU to divide the grid integration into subtasks.
 * 2. Copy the subtask information to the GPU.
 * 3. Calculate the matrix elements on the GPU.
 * 4. Perform matrix multiplication on the GPU.
 * 5. Copy the results back to the host.
 */
void gint_vl_gpu(hamilt::HContainer<double>* hRGint,
                 const double* vlocal,
                 const double* ylmcoef_now,
                 const double dr,
                 const double* rcut,
                 const Grid_Technique& gridt,
                 const UnitCell& ucell)
{
    checkCuda(cudaSetDevice(gridt.dev_id));
    // checkCuda(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    const int nbzp = gridt.nbzp;
    const int num_streams = gridt.nstreams;
    const int max_atom = gridt.max_atom;
    const int max_atom_per_bcell = max_atom * gridt.bxyz;
    const int max_atom_per_z = max_atom_per_bcell * nbzp;
    const int max_phi_per_z = max_atom_per_z * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;
    const double vfactor = ucell.omega / gridt.ncxyz;
    const int nczp = nbzp * gridt.bz;
    std::vector<cudaStream_t> streams(num_streams);

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }

    const int nnrg = hRGint->get_nnr();
    hRGint->set_zero();
    Cuda_Mem_Wrapper<double> grid_vlocal_g(nnrg, 1, false);
    grid_vlocal_g.memset_device_sync();

    Cuda_Mem_Wrapper<double> dr_part(max_atom_per_z * 3, num_streams, true);
    Cuda_Mem_Wrapper<uint8_t> atoms_type(max_atom_per_z, num_streams, true);
    // The first number in every group of two represents the number of atoms on that bigcell.
    // The second number represents the cumulative number of atoms up to that bigcell.
    Cuda_Mem_Wrapper<int> atoms_num_info(2 * nbzp, num_streams, true);
    Cuda_Mem_Wrapper<double> vldr3(nbzp * gridt.bxyz, num_streams, true);

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

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_streams) collapse(2)
#endif
    for (int i = 0; i < gridt.nbx; i++)
    {
        for (int j = 0; j < gridt.nby; j++)
        {
            // 20240620 Note that it must be set again here because 
            // cuda's device is not safe in a multi-threaded environment.
            checkCuda(cudaSetDevice(gridt.dev_id));
#ifdef _OPENMP
            const int sid = omp_get_thread_num();
#else
            const int sid = 0;
#endif

            int max_m = 0;
            int max_n = 0;
            int atom_pair_num = 0;
            int atoms_per_z = 0;
            const int grid_index_ij = i * gridt.nby * nbzp + j * nbzp;
            
            gtask_vlocal(gridt,
                         ucell,
                         grid_index_ij,
                         nczp,
                         vfactor,
                         vlocal,
                         atoms_per_z,
                         atoms_num_info.get_host_pointer(sid),
                         atoms_type.get_host_pointer(sid),
                         dr_part.get_host_pointer(sid),
                         vldr3.get_host_pointer(sid));
        
            alloc_mult_vlocal(hRGint,
                              gridt,
                              ucell,
                              grid_index_ij,
                              max_atom,
                              psi.get_device_pointer(sid),
                              psi_vldr3.get_device_pointer(sid),
                              grid_vlocal_g.get_device_pointer(),
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

            dr_part.copy_host_to_device_async(streams[sid], sid, atoms_per_z * 3);
            atoms_type.copy_host_to_device_async(streams[sid], sid, atoms_per_z);
            vldr3.copy_host_to_device_async(streams[sid], sid);
            atoms_num_info.copy_host_to_device_async(streams[sid], sid, 2 * nbzp);
            
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
            psi_vldr3.memset_device_async(streams[sid], sid, 0);

            dim3 grid_psi(nbzp, gridt.bxyz);
            dim3 block_psi(64);
            get_psi_and_vldr3<<<grid_psi,
                                block_psi,
                                0,
                                streams[sid]>>>(
                gridt.ylmcoef_g,
                dr,
                gridt.bxyz,
                ucell.nwmax,
                max_atom,
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_nw_g,
                gridt.rcut_g,
                gridt.nr_max,
                gridt.psi_u_g,
                gridt.mcell_pos_g,
                dr_part.get_device_pointer(sid),
                vldr3.get_device_pointer(sid),
                atoms_type.get_device_pointer(sid),
                atoms_num_info.get_device_pointer(sid),
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
            checkCuda(cudaStreamSynchronize(streams[sid]));
        }
    }

    checkCuda(cudaMemcpy(
        hRGint->get_wrapper(),
        grid_vlocal_g.get_device_pointer(),
        nnrg * sizeof(double),
        cudaMemcpyDeviceToHost));

    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamDestroy(streams[i]));
    }
}

} // namespace GintKernel
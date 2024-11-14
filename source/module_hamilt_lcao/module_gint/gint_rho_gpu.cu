#include "kernels/cuda/cuda_tools.cuh"
#include "module_base/ylm.h"
#include "gint_rho_gpu.h"
#include "gint_tools.h"
#include "kernels/cuda/gint_rho.cuh"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GintKernel
{

void gint_rho_gpu(const hamilt::HContainer<double>* dm,
                        const double* ylmcoef_now,
                        const double dr,
                        const double* rcut,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        double* rho)
{
    checkCuda(cudaSetDevice(gridt.dev_id));
    // checkCuda(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));

    const int nbzp = gridt.nbzp;
    const int nczp =nbzp * gridt.bz;
    const int num_mcell_on_proc = nczp * gridt.ncx * gridt.ncy;
    const int lgd = gridt.lgd;
    const int max_atom = gridt.max_atom;
    const int num_streams = gridt.nstreams;
    const int max_atom_per_bcell = max_atom * gridt.bxyz;
    const int max_atom_per_z = max_atom * nbzp;
    const int max_phi_per_z = max_atom_per_bcell * nbzp * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;

    std::vector<cudaStream_t> streams(num_streams);
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }

    Cuda_Mem_Wrapper<double> dr_part(max_atom_per_z * 3, num_streams, true);
    Cuda_Mem_Wrapper<uint8_t> atoms_type(max_atom_per_z, num_streams, true);
    // The first number in every group of two represents the number of atoms on that bigcell.
    // The second number represents the cumulative number of atoms up to that bigcell.
    Cuda_Mem_Wrapper<int> atoms_num_info(2 * nbzp, num_streams, true);

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

    Cuda_Mem_Wrapper<double> dm_matrix(dm->get_nnr(), 1, false);
    // retrieve the density matrix on the host
    checkCuda(cudaMemcpy(dm_matrix.get_device_pointer(),
                         dm->get_wrapper(),
                         dm->get_nnr() * sizeof(double),
                         cudaMemcpyHostToDevice));

// calculate the rho for every nbzp bigcells
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
            // get stream id
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

            // generate GPU tasks, including the calculation of psir, matrix
            // multiplication, and dot product
            gtask_rho(gridt,
                      grid_index_ij,
                      ucell,
                      dr_part.get_host_pointer(sid),
                      atoms_type.get_host_pointer(sid),
                      atoms_num_info.get_host_pointer(sid),
                      atoms_per_z);

            alloc_mult_dot_rho(
                dm,
                gridt,
                ucell,
                grid_index_ij,
                max_atom,
                lgd,
                nczp,
                atoms_num_info.get_host_pointer(sid),
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
            
            dr_part.copy_host_to_device_async(streams[sid], sid, atoms_per_z * 3);
            atoms_type.copy_host_to_device_async(streams[sid], sid, atoms_per_z);
            atoms_num_info.copy_host_to_device_async(streams[sid], sid);

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
            dim3 grid_psi(nbzp, gridt.bxyz);
            dim3 block_psi(64);
            get_psi<<<grid_psi, block_psi, 0, streams[sid]>>>(
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
                atoms_type.get_device_pointer(sid),
                atoms_num_info.get_device_pointer(sid),
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
            // if warpSize is not eauql to 32, the psir_dot kernel should be modified
            dim3 grid_dot(nbzp, gridt.bxyz);
            dim3 block_dot(64); 
            psir_dot<<<grid_dot, block_dot, sizeof(double) * 32, streams[sid]>>>(
                gridt.bxyz,
                ucell.nwmax,
                atoms_num_info.get_device_pointer(sid),
                psi.get_device_pointer(sid),
                psi_dm.get_device_pointer(sid),
                dot_product.get_device_pointer(sid));
            checkCudaLastError();
            checkCuda(cudaStreamSynchronize(streams[sid]));
        }
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

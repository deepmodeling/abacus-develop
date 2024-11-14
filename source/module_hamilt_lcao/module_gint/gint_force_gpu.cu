#ifdef _OPENMP
#include <omp.h>
#endif

#include "gint_force_gpu.h"
#include "kernels/cuda/cuda_tools.cuh"
#include "kernels/cuda/gint_force.cuh"
#include "module_base/ylm.h"
#include "gint_tools.h"

namespace GintKernel
{
/**
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
void gint_fvl_gpu(const hamilt::HContainer<double>* dm,
                        const double* vlocal,
                        double* force_in,
                        double* stress_in,
                        double dr,
                        const double* rcut,
                        const int isforce,
                        const int isstress,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell)
{ 
    checkCuda(cudaSetDevice(gridt.dev_id));
    // checkCuda(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));

    const int nbzp = gridt.nbzp;
    const int max_atom = gridt.max_atom;
    const int nwmax = ucell.nwmax;
    const int bxyz = gridt.bxyz;
    const int max_atom_per_bcell = max_atom * bxyz;
    const int max_atom_per_z = max_atom_per_bcell * nbzp;
    const int max_phi_per_z = max_atom_per_z * ucell.nwmax;
    const int max_atompair_per_z = max_atom * max_atom * nbzp;
    const double vfactor = ucell.omega / gridt.ncxyz;
    const int nczp = nbzp * gridt.bz;
    const int nat=ucell.nat;

    const int num_streams = gridt.nstreams;

    std::vector<cudaStream_t> streams(num_streams);
    for (int i = 0; i < num_streams; i++)
    {
        checkCuda(cudaStreamCreate(&streams[i]));
    }

    Cuda_Mem_Wrapper<double> dr_part(3 * max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<uint8_t> atoms_type(max_atom_per_z, num_streams, true);
    Cuda_Mem_Wrapper<int> iat_on_nbz(max_atom_per_z, num_streams, true);
    // The first number in every group of two represents the number of atoms on that bigcell.
    // The second number represents the cumulative number of atoms up to that bigcell.
    Cuda_Mem_Wrapper<int> atoms_num_info(2 * nbzp, num_streams, true);
    Cuda_Mem_Wrapper<double> vldr3(nbzp * gridt.bxyz, num_streams, true);

    Cuda_Mem_Wrapper<double> psi(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> psi_dm(max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> dpsi(3 * max_phi_per_z, num_streams, false);
    Cuda_Mem_Wrapper<double> d2psi(6 * max_phi_per_z, num_streams, false);

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

    Cuda_Mem_Wrapper<double> dm_matrix(dm->get_nnr(), 1, false);
    // retrieve the density matrix on the host
    checkCuda(cudaMemcpy(dm_matrix.get_device_pointer(),
                         dm->get_wrapper(),
                         dm->get_nnr() * sizeof(double),
                         cudaMemcpyHostToDevice));
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

            gtask_force(gridt,
                        ucell,
                        grid_index_ij,
                        nczp,
                        vfactor,
                        vlocal,
                        atoms_per_z,
                        atoms_num_info.get_host_pointer(sid),
                        iat_on_nbz.get_host_pointer(sid),
                        atoms_type.get_host_pointer(sid),
                        dr_part.get_host_pointer(sid),
                        vldr3.get_host_pointer(sid));
           
            alloc_mult_force(dm,
                             gridt,
                             ucell, 
                             grid_index_ij,
                             max_atom,
                             atoms_num_info.get_host_pointer(sid),
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
                             gemm_C.get_host_pointer(sid));

            dr_part.copy_host_to_device_async(streams[sid], sid, 3 * atoms_per_z);
            atoms_type.copy_host_to_device_async(streams[sid], sid, atoms_per_z);
            iat_on_nbz.copy_host_to_device_async(streams[sid], sid, atoms_per_z);
            vldr3.copy_host_to_device_async(streams[sid], sid);
            atoms_num_info.copy_host_to_device_async(streams[sid], sid);
            
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
            dpsi.memset_device_async(streams[sid], sid, 0);
            d2psi.memset_device_async(streams[sid], sid, 0);

            dim3 grid_psi(nbzp, gridt.bxyz);
            dim3 block_psi(64);
            get_psi_force<<<grid_psi,
                            block_psi,
                            0,
                            streams[sid]>>>(
                gridt.ylmcoef_g,
                dr,
                bxyz,
                nwmax,
                max_atom,
                gridt.atom_nwl_g,
                gridt.atom_new_g,
                gridt.atom_ylm_g,
                gridt.atom_l_g,
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
                dpsi.get_device_pointer(sid),
                d2psi.get_device_pointer(sid));
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
                dim3 grid_force(nbzp);
                dim3 block_force(64);
                dot_product_force<<<grid_force,
                                    block_force,
                                    32 * 3 * sizeof(double),
                                    streams[sid]>>>(
                                        bxyz,
                                        nwmax,
                                        atoms_num_info.get_device_pointer(sid),
                                        iat_on_nbz.get_device_pointer(sid),
                                        dpsi.get_device_pointer(sid),
                                        psi_dm.get_device_pointer(sid),
                                        force.get_device_pointer(sid));
                checkCudaLastError();
            }

            if (isstress){ 
                dim3 grid_stress(nbzp);
                dim3 block_stress(64);
                dot_product_stress<<<grid_stress,
                                        block_stress,
                                        32 * 6 * sizeof(double),
                                        streams[sid]>>>(
                                    d2psi.get_device_pointer(sid),
                                    psi_dm.get_device_pointer(sid),
                                    atoms_per_z * nwmax * bxyz,
                                    stress.get_device_pointer(sid));
                checkCudaLastError();
            }
            checkCuda(cudaStreamSynchronize(streams[sid]));
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

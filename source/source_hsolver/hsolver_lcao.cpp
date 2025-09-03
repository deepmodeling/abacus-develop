#include "hsolver_lcao.h"

#ifdef __MPI
#include "diago_scalapack.h"
#include "source_base/module_external/scalapack_connector.h"
#else
#include "diago_lapack.h"
#endif

#ifdef __CUSOLVERMP
#include "diago_cusolvermp.h"
#endif

#ifdef __ELPA
#include "diago_elpa.h"
#include "diago_elpa_native.h"
#endif

#ifdef __CUDA
#include "diago_cusolver.h"
#endif

#ifdef __PEXSI
#include "diago_pexsi.h"
#endif

#include "source_base/global_variable.h"
#include "source_estate/elecstate_tools.h"
#include "source_base/memory.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hsolver/parallel_k2d.h"
#include "source_io/module_parameter/parameter.h"

namespace hsolver
{

template <typename T, typename Device>
void HSolverLCAO<T, Device>::solve(hamilt::Hamilt<T>* pHamilt,
                                   psi::Psi<T>& psi,
                                   elecstate::ElecState* pes,
                                   const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLCAO", "solve");
    ModuleBase::timer::tick("HSolverLCAO", "solve");

    if (this->method != "pexsi")
    {
    #ifdef __MPI
    #ifdef __CUDA
        if (this->method == "cusolver" && GlobalV::NPROC > 1)
        {
            this->parakSolve_cusolver(pHamilt, psi, pes);
        }else 
    #endif
        if (PARAM.globalv.kpar_lcao > 1
            && (this->method == "genelpa" || this->method == "elpa" || this->method == "scalapack_gvx"))
        {
            this->parakSolve(pHamilt, psi, pes, PARAM.globalv.kpar_lcao);
        } else
    #endif
        if (PARAM.globalv.kpar_lcao == 1)
        {
            /// Loop over k points for solve Hamiltonian to eigenpairs(eigenvalues and eigenvectors).
            for (int ik = 0; ik < psi.get_nk(); ++ik)
            {
                /// update H(k) for each k point
                pHamilt->updateHk(ik);

                /// find psi pointer for each k point
                psi.fix_k(ik);

                /// solve eigenvector and eigenvalue for H(k)
                this->hamiltSolvePsiK(pHamilt, psi, &(pes->ekb(ik, 0)));
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("HSolverLCAO::solve",
                                     "This method and KPAR setting is not supported for lcao basis in ABACUS!");
        }

        elecstate::calculate_weights(pes->ekb,
                                     pes->wg,
                                     pes->klist,
                                     pes->eferm,
                                     pes->f_en,
                                     pes->nelec_spin,
                                     pes->skip_weights);

        auto _pes_lcao = dynamic_cast<elecstate::ElecStateLCAO<T>*>(pes);
        elecstate::calEBand(_pes_lcao->ekb, _pes_lcao->wg, _pes_lcao->f_en);
        elecstate::cal_dm_psi(_pes_lcao->DM->get_paraV_pointer(), _pes_lcao->wg, psi, *(_pes_lcao->DM));
        _pes_lcao->DM->cal_DMR();

        if (!skip_charge)
        {
            // used in scf calculation
            // calculate charge by eigenpairs(eigenvalues and eigenvectors)
            pes->psiToRho(psi);
        }
        else
        {
            // used in nscf calculation
        }
    }
    else if (this->method == "pexsi")
    {
#ifdef __PEXSI // other purification methods should follow this routine
        DiagoPexsi<T> pe(ParaV);
        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            /// update H(k) for each k point
            pHamilt->updateHk(ik);
            psi.fix_k(ik);
            // solve eigenvector and eigenvalue for H(k)
            pe.diag(pHamilt, psi, nullptr);
        }
        auto _pes = dynamic_cast<elecstate::ElecStateLCAO<T>*>(pes);
        pes->f_en.eband = pe.totalFreeEnergy;
        // maybe eferm could be dealt with in the future
        _pes->dmToRho(pe.DM, pe.EDM);
#endif
    }

    ModuleBase::timer::tick("HSolverLCAO", "solve");
    return;
}

template <typename T, typename Device>
void HSolverLCAO<T, Device>::hamiltSolvePsiK(hamilt::Hamilt<T>* hm, psi::Psi<T>& psi, double* eigenvalue)
{
    ModuleBase::TITLE("HSolverLCAO", "hamiltSolvePsiK");
    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");

    if (this->method == "scalapack_gvx")
    {
#ifdef __MPI
        DiagoScalapack<T> sa;
        sa.diag(hm, psi, eigenvalue);
#endif
    }
#ifdef __ELPA
    else if (this->method == "genelpa")
    {
        DiagoElpa<T> el;
        el.diag(hm, psi, eigenvalue);
    }
    else if (this->method == "elpa")
    {
        DiagoElpaNative<T> el;
        el.diag(hm, psi, eigenvalue);
    }
#endif
#ifdef __CUDA
    else if (this->method == "cusolver")
    {
        // Note: This branch will only be executed in the single-process case
        DiagoCusolver<T> cu;
        hamilt::MatrixBlock<T> hk, sk;
        hm->matrix(hk, sk);
        cu.diag(hk, sk, psi, eigenvalue);
    }
#ifdef __CUSOLVERMP
    else if (this->method == "cusolvermp")
    {
        DiagoCusolverMP<T> cm;
        cm.diag(hm, psi, eigenvalue);
    }
#endif
#endif
#ifndef __MPI
    else if (this->method == "lapack") // only for single core
    {
        DiagoLapack<T> la;
        la.diag(hm, psi, eigenvalue);
    }
#endif
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method is not supported for lcao basis in ABACUS!");
    }

    ModuleBase::timer::tick("HSolverLCAO", "hamiltSolvePsiK");
}

template <typename T, typename Device>
void HSolverLCAO<T, Device>::parakSolve(hamilt::Hamilt<T>* pHamilt,
                                        psi::Psi<T>& psi,
                                        elecstate::ElecState* pes,
                                        int kpar)
{
#ifdef __MPI
    ModuleBase::timer::tick("HSolverLCAO", "parakSolve");
    auto k2d = Parallel_K2D<T>();
    k2d.set_kpar(kpar);
    int nbands = this->ParaV->get_nbands();
    int nks = psi.get_nk();
    int nrow = this->ParaV->get_global_row_size();
    int nb2d = this->ParaV->get_block_size();
    k2d.set_para_env(psi.get_nk(), nrow, nb2d, GlobalV::NPROC, GlobalV::MY_RANK, PARAM.inp.nspin);
    /// set psi_pool
    const int zero = 0;
    int ncol_bands_pool
        = numroc_(&(nbands), &(nb2d), &(k2d.get_p2D_pool()->coord[1]), &zero, &(k2d.get_p2D_pool()->dim1));
    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < k2d.get_pKpoints()->get_max_nks_pool(); ++ik)
    {
        // if nks is not equal to the number of k points in the pool
        std::vector<int> ik_kpar;
        int ik_avail = 0;
        for (int i = 0; i < k2d.get_kpar(); i++)
        {
            if (ik + k2d.get_pKpoints()->startk_pool[i] < nks && ik < k2d.get_pKpoints()->nks_pool[i])
            {
                ik_avail++;
            }
        }
        if (ik_avail == 0)
        {
            ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "ik_avail is 0!");
        }
        else
        {
            ik_kpar.resize(ik_avail);
            for (int i = 0; i < ik_avail; i++)
            {
                ik_kpar[i] = ik + k2d.get_pKpoints()->startk_pool[i];
            }
        }
        k2d.distribute_hsk(pHamilt, ik_kpar, nrow);
        /// global index of k point
        int ik_global = ik + k2d.get_pKpoints()->startk_pool[k2d.get_my_pool()];
        auto psi_pool = psi::Psi<T>(1, ncol_bands_pool, k2d.get_p2D_pool()->nrow, k2d.get_p2D_pool()->nrow, true);
        ModuleBase::Memory::record("HSolverLCAO::psi_pool", nrow * ncol_bands_pool * sizeof(T));
        if (ik_global < psi.get_nk() && ik < k2d.get_pKpoints()->nks_pool[k2d.get_my_pool()])
        {
            /// local psi in pool
            psi_pool.fix_k(0);
            hamilt::MatrixBlock<T> hk_pool = hamilt::MatrixBlock<T>{k2d.hk_pool.data(),
                                                                    (size_t)k2d.get_p2D_pool()->get_row_size(),
                                                                    (size_t)k2d.get_p2D_pool()->get_col_size(),
                                                                    k2d.get_p2D_pool()->desc};
            hamilt::MatrixBlock<T> sk_pool = hamilt::MatrixBlock<T>{k2d.sk_pool.data(),
                                                                    (size_t)k2d.get_p2D_pool()->get_row_size(),
                                                                    (size_t)k2d.get_p2D_pool()->get_col_size(),
                                                                    k2d.get_p2D_pool()->desc};
            /// solve eigenvector and eigenvalue for H(k)
            if (this->method == "scalapack_gvx")
            {
                DiagoScalapack<T> sa;
                sa.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
#ifdef __ELPA
            else if (this->method == "genelpa")
            {
                DiagoElpa<T> el;
                el.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
            else if (this->method == "elpa")
            {
                DiagoElpaNative<T> el;
                el.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
#endif
            else
            {
                ModuleBase::WARNING_QUIT("HSolverLCAO::solve",
                                         "This type of eigensolver for k-parallelism diagnolization is not supported!");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ModuleBase::timer::tick("HSolverLCAO", "collect_psi");
        for (int ipool = 0; ipool < ik_kpar.size(); ++ipool)
        {
            int source = k2d.get_pKpoints()->get_startpro_pool(ipool);
            MPI_Bcast(&(pes->ekb(ik_kpar[ipool], 0)), nbands, MPI_DOUBLE, source, MPI_COMM_WORLD);
            int desc_pool[9];
            std::copy(k2d.get_p2D_pool()->desc, k2d.get_p2D_pool()->desc + 9, desc_pool);
            if (k2d.get_my_pool() != ipool)
            {
                desc_pool[1] = -1;
            }
            psi.fix_k(ik_kpar[ipool]);
            Cpxgemr2d(nrow,
                      nbands,
                      psi_pool.get_pointer(),
                      1,
                      1,
                      desc_pool,
                      psi.get_pointer(),
                      1,
                      1,
                      k2d.get_p2D_global()->desc,
                      k2d.get_p2D_global()->blacs_ctxt);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ModuleBase::timer::tick("HSolverLCAO", "collect_psi");
    }
    k2d.unset_para_env();
    ModuleBase::timer::tick("HSolverLCAO", "parakSolve");
#endif
}

#if defined (__MPI) && defined (__CUDA)
template <typename T, typename Device>
void HSolverLCAO<T, Device>::parakSolve_cusolver(hamilt::Hamilt<T>* pHamilt,
                                            psi::Psi<T>& psi,
                                            elecstate::ElecState* pes)
{
    ModuleBase::timer::tick("HSolverLCAO", "parakSolve");
    const int dev_id = base_device::information::set_device_by_rank();

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Split communicator by shared memory node
    MPI_Comm nodeComm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &nodeComm);

    int local_rank, local_size;
    MPI_Comm_rank(nodeComm, &local_rank);
    MPI_Comm_size(nodeComm, &local_size);

    // Get number of CUDA devices on this node
    int device_count = 0;
    cudaError_t cuda_err = cudaGetDeviceCount(&device_count);
    if (cuda_err != cudaSuccess) {
        device_count = 0; // Treat as no GPU available
    }

    if(local_rank >= device_count) {
        local_rank = -1; // Mark as inactive for GPU work
    }

    // Determine the number of MPI processes on this node that can actively use a GPU.
    // This is the minimum of:
    //   - The number of available MPI processes on the node (local_size)
    //   - The number of available CUDA-capable GPUs on the node (device_count)
    // Each GPU is assumed to be used by one dedicated MPI process.
    // Thus, only the first 'min(local_size, device_count)' ranks on this node
    // will be assigned GPU work; the rest will be inactive or used for communication-only roles.
    int active_procs_per_node = std::min(local_size, device_count);

    std::vector<int> all_active_procs(world_size);
    std::vector<int> all_local_ranks(world_size);

    MPI_Allgather(&active_procs_per_node, 1, MPI_INT, 
                  all_active_procs.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&local_rank, 1, MPI_INT, 
                  all_local_ranks.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    int total_active_ranks = 0;
    int total_nodes = 0;
    int highest_active_rank = 0;

    for (int i = 0; i < world_size; ++i) {
        if (all_local_ranks[i] == 0) {  // new node
            total_nodes++;
            total_active_ranks += all_active_procs[i];
            highest_active_rank = std::max(highest_active_rank, all_active_procs[i] - 1);
        }
    }

    // active_ranks will store the global ranks of all active processes across all nodes
    // The order of global ranks stored here determines the order in which they will be assigned to k-points.
    // The k-points will be distributed among these ranks in a round-robin fashion.
    // The purpose of setting the order is to ensure load balancing among nodes as much as possible
    std::vector<int> active_ranks;
    for(int i = 0; i <= highest_active_rank; i++)
    {
        for(int j = 0; j < world_size; j++)
        {
            if(all_local_ranks[j] == i)
            {
                active_ranks.push_back(j);
            }
        }
    }
    
    const int nks = psi.get_nk();  // total number of k points
    const int nbands = this->ParaV->get_nbands();
    // Set the parallel storage scheme for the matrix and psi
    Parallel_2D mat_para_global;    // store the info about how the origin matrix is distributed in parallel
    Parallel_2D mat_para_local;     // store the info about how the matrix is distributed after collected from all processes
    Parallel_2D psi_para_global;    // store the info about how the psi is distributed in parallel
    Parallel_2D psi_para_local;     // store the info about how the psi is distributed before distributing to all processes

    MPI_Comm self_comm;  // the communicator that only contains the current process itself
    MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &self_comm);
    int nrow = this->ParaV->get_global_row_size(); // number of rows in the global matrix
    int ncol = nrow;
    int nb2d = this->ParaV->get_block_size();      // block size for the 2D matrix distribution
    mat_para_global.init(nrow, ncol, nb2d, MPI_COMM_WORLD);
    psi_para_global.init(nrow, nbands, nb2d, MPI_COMM_WORLD);
    mat_para_local.init(nrow, ncol, nb2d, self_comm);
    psi_para_local.init(nrow, ncol, nb2d, self_comm);

    std::vector<T> hk_mat; // temporary storage for H(k) matrix collected from all processes
    std::vector<T> sk_mat; // temporary storage for S(k) matrix collected from all processes
    // In each iteration, we process total_active_ranks k-points.
    for(int ik_start = 0; ik_start < nks; ik_start += total_active_ranks)
    {
        int kpt_assigned = -1; // the k-point assigned to the current MPI process in this iteration
        // Compute and gather the hk and sk matrices distributed across different processes in parallel,
        // preparing for subsequent transfer to the GPU for computation.
        for(int ik = ik_start; ik < ik_start + total_active_ranks && ik < nks; ik++)
        {
            // `is_active` indicates whether this MPI process is assigned to compute the current k-point
            bool is_active = world_rank == active_ranks[ik % total_active_ranks];
            if (is_active)
            {
                kpt_assigned = ik;
                hk_mat.resize(nrow * ncol);
                sk_mat.resize(nrow * ncol);
            }
            pHamilt->updateHk(ik);
            hamilt::MatrixBlock<T> hk_2D, sk_2D;
            pHamilt->matrix(hk_2D, sk_2D);
            int desc_tmp[9];
            T* hk_local_ptr = hk_mat.data();
            T* sk_local_ptr = sk_mat.data();
            std::copy(mat_para_local.desc, mat_para_local.desc + 9, desc_tmp);
            if( !is_active)
            {
                desc_tmp[1] = -1;
            }

            Cpxgemr2d(nrow, ncol, hk_2D.p, 1, 1, mat_para_global.desc,
                      hk_local_ptr, 1, 1, desc_tmp,
                      mat_para_global.blacs_ctxt);
            Cpxgemr2d(nrow, ncol, sk_2D.p, 1, 1, mat_para_global.desc,
                      sk_local_ptr, 1, 1, desc_tmp,
                      mat_para_global.blacs_ctxt);
        }

        // diagonalize the Hamiltonian matrix using cusolver
        psi::Psi<T> psi_local{};
        if(kpt_assigned != -1)
        {
            psi_local.resize(1, ncol, nrow);
            DiagoCusolver<T> cu{};
            hamilt::MatrixBlock<T> hk_local = hamilt::MatrixBlock<T>{
                    hk_mat.data(), (size_t)nrow, (size_t)ncol,
                    mat_para_local.desc};
            hamilt::MatrixBlock<T> sk_local = hamilt::MatrixBlock<T>{
                    sk_mat.data(), (size_t)nrow, (size_t)ncol,
                    mat_para_local.desc};
            cu.diag(hk_local, sk_local, psi_local, &(pes->ekb(kpt_assigned, 0)));
        }

        // transfer the eigenvectors and eigenvalues to all processes
        for(int ik = ik_start; ik < ik_start + total_active_ranks && ik < nks; ik++)
        {
            int root = active_ranks[ik % total_active_ranks];
            MPI_Bcast(&(pes->ekb(ik, 0)), nbands, MPI_DOUBLE, root, MPI_COMM_WORLD);
            int desc_pool[9];
            std::copy(psi_para_local.desc, psi_para_local.desc + 9, desc_pool);
            T* psi_local_ptr = nullptr;
            if (world_rank != root)
            {
                desc_pool[1] = -1;
            }else
            {
                psi_local_ptr = psi_local.get_pointer();
            }
            psi.fix_k(ik);
            Cpxgemr2d(nrow,
                      nbands,
                      psi_local_ptr,
                      1,
                      1,
                      desc_pool,
                      psi.get_pointer(),
                      1,
                      1,
                      psi_para_global.desc,
                      psi_para_global.blacs_ctxt);
        }
    }

    MPI_Comm_free(&self_comm);
    MPI_Comm_free(&nodeComm);
    ModuleBase::timer::tick("HSolverLCAO", "parakSolve");
}
#endif


template class HSolverLCAO<double>;
template class HSolverLCAO<std::complex<double>>;

} // namespace hsolver
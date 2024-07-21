#include "parallel_k2d.h"

#include "module_base/parallel_global.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_base/memory.h"

template <typename TK>
void Parallel_K2D<TK>::set_para_env(int nks,
                                    const int& nw,
                                    const int& nb2d,
                                    const int& nproc,
                                    const int& my_rank,
                                    const int& nspin) {
    const int kpar = this->get_kpar();
    Parallel_Global::divide_mpi_groups(nproc,
                                       kpar,
                                       my_rank,
                                       this->NPROC_IN_POOL,
                                       this->MY_POOL,
                                       this->RANK_IN_POOL);
#ifdef __MPI
    MPI_Comm_split(MPI_COMM_WORLD,
                   this->MY_POOL,
                   this->RANK_IN_POOL,
                   &this->POOL_WORLD_K2D);
#endif
    this->Pkpoints = new Parallel_Kpoints;
    this->P2D_global = new Parallel_2D;
    this->P2D_pool = new Parallel_2D;
    this->Pkpoints
        ->kinfo(nks, kpar, this->MY_POOL, this->RANK_IN_POOL, nproc, nspin);
    this->P2D_global->init(nw, nw, nb2d, MPI_COMM_WORLD);
    this->P2D_pool->init(nw, nw, nb2d, this->POOL_WORLD_K2D);
}

template <typename TK>
void Parallel_K2D<TK>::distribute_hsk(hamilt::Hamilt<TK>* pHamilt,
                                      const std::vector<int>& ik_kpar,
                                      const int& nw) {
#ifdef __MPI
    ModuleBase::timer::tick("Parallel_K2D", "distribute_hsk");
    for (int ipool = 0; ipool < ik_kpar.size(); ++ipool)
    {
        pHamilt->updateHk(ik_kpar[ipool]);
        hamilt::MatrixBlock<TK> HK_global, SK_global;
        pHamilt->matrix(HK_global, SK_global);
        if (this->MY_POOL == this->Pkpoints->whichpool[ik_kpar[ipool]]) {
            this->hk_pool.resize(this->P2D_pool->get_local_size(), 0.0);
            this->sk_pool.resize(this->P2D_pool->get_local_size(), 0.0);
        }
        int desc_pool[9];
        std::copy(this->P2D_pool->desc, this->P2D_pool->desc + 9, desc_pool);
        if (this->MY_POOL != this->Pkpoints->whichpool[ik_kpar[ipool]]) {
            desc_pool[1] = -1;
        }
        Cpxgemr2d(nw,
                  nw,
                  HK_global.p,
                  1,
                  1,
                  this->P2D_global->desc,
                  hk_pool.data(),
                  1,
                  1,
                  desc_pool,
                  this->P2D_global->blacs_ctxt);
        Cpxgemr2d(nw,
                  nw,
                  SK_global.p,
                  1,
                  1,
                  this->P2D_global->desc,
                  sk_pool.data(),
                  1,
                  1,
                  desc_pool,
                  this->P2D_global->blacs_ctxt);
    }
    ModuleBase::Memory::record("Parallel_K2D::hsk_pool", this->P2D_pool->get_local_size() * 2 * sizeof(TK));
    ModuleBase::timer::tick("Parallel_K2D", "distribute_hsk");
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

template <typename TK>
void Parallel_K2D<TK>::unset_para_env() {
    if (this->Pkpoints != nullptr) {
        delete this->Pkpoints;
        this->Pkpoints = nullptr;
    }
    if (this->P2D_global != nullptr) {
        delete this->P2D_global;
        this->P2D_global = nullptr;
    }
    if (this->P2D_pool != nullptr) {
        delete this->P2D_pool;
        this->P2D_pool = nullptr;
    }
    MPI_Comm_free(&this->POOL_WORLD_K2D);
}

template <typename TK>
void Parallel_K2D<TK>::set_kpar(int kpar) {
    if (kpar < 1) {
        ModuleBase::WARNING_QUIT("Parallel_K2D::set_kpar",
                                 "kpar must be greater than 0.");
    }
    this->kpar_ = kpar;
}

template class Parallel_K2D<double>;
template class Parallel_K2D<std::complex<double>>;
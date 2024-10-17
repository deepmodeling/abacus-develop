#include "hsolver_lcao.h"

#ifdef __MPI
#include "diago_scalapack.h"
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
#include "module_elecstate/elecstate_lcao.h"
#endif

#include "diago_cg.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hsolver/parallel_k2d.h"
#include "module_io/write_HS.h"
#include "module_parameter/parameter.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#include <unistd.h>

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

#ifdef __PEXSI // other purification methods should follow this routine
    // Zhang Xiaoyang :  Please modify Pesxi usage later
    if (this->method == "pexsi")
    {
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
        ModuleBase::timer::tick("HSolverLCAO", "solve");
        return;
    }
#endif

    if (this->method == "cg_in_lcao")
    {
        this->precondition_lcao.resize(psi.get_nbasis());

        using Real = typename GetTypeReal<T>::type;
        // set precondition
        for (size_t i = 0; i < precondition_lcao.size(); i++)
        {
            precondition_lcao[i] = 1.0;
        }
    }

#ifdef __MPI
    if (GlobalV::KPAR_LCAO > 1
        && (this->method == "genelpa" || this->method == "elpa" || this->method == "scalapack_gvx"))
    {
        this->parakSolve(pHamilt, psi, pes, GlobalV::KPAR_LCAO);
    }
    else
#endif
    {
        /// Loop over k points for solve Hamiltonian to charge density
        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            /// update H(k) for each k point
            pHamilt->updateHk(ik);

            psi.fix_k(ik);

            // solve eigenvector and eigenvalue for H(k)
            this->hamiltSolvePsiK(pHamilt, psi, &(pes->ekb(ik, 0)));
        }
    }

    if (skip_charge) // used in nscf calculation
    {
        ModuleBase::timer::tick("HSolverLCAO", "solve");
    }
    else // used in scf calculation
    {
        // calculate charge by psi
        pes->psiToRho(psi);
        ModuleBase::timer::tick("HSolverLCAO", "solve");
    }

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
        DiagoCusolver<T> cs(this->ParaV);
        cs.diag(hm, psi, eigenvalue);
    }
    else if (this->method == "cusolvermp")
    {
#ifdef __CUSOLVERMP
        DiagoCusolverMP<T> cm;
        cm.diag(hm, psi, eigenvalue);
#else
        ModuleBase::WARNING_QUIT("HSolverLCAO", "CUSOLVERMP did not compiled!");
#endif
    }
#endif
    else if (this->method == "lapack")
    {
#ifndef __MPI
        DiagoLapack<T> la;
        la.diag(hm, psi, eigenvalue);
#else
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This type of eigensolver is not supported!");
#endif
    }
    else
    {

        using ct_Device = typename ct::PsiToContainer<base_device::DEVICE_CPU>::type;

        auto subspace_func = [](const ct::Tensor& psi_in, ct::Tensor& psi_out) {
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim == 2, "dims of psi_in should be less than or equal to 2");
        };

        DiagoCG<T, Device> cg(PARAM.inp.basis_type,
                              PARAM.inp.calculation,
                              DiagoIterAssist<T, Device>::need_subspace,
                              subspace_func,
                              DiagoIterAssist<T, Device>::PW_DIAG_THR,
                              DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                              GlobalV::NPROC_IN_POOL);

        hamilt::MatrixBlock<T> h_mat, s_mat;
        hm->matrix(h_mat, s_mat);

        // set h_mat & s_mat
        for (int i = 0; i < h_mat.row; i++)
        {
            for (int j = i; j < h_mat.col; j++)
            {
                h_mat.p[h_mat.row * j + i] = hsolver::get_conj(h_mat.p[h_mat.row * i + j]);
                s_mat.p[s_mat.row * j + i] = hsolver::get_conj(s_mat.p[s_mat.row * i + j]);
            }
        }

        const T *one_ = nullptr, *zero_ = nullptr;
        one_ = new T(static_cast<T>(1.0));
        zero_ = new T(static_cast<T>(0.0));

        auto hpsi_func = [h_mat, one_, zero_](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");

            Device* ctx = {};

            gemv_op<T, Device>()(ctx,
                                 'N',
                                 h_mat.row,
                                 h_mat.col,
                                 one_,
                                 h_mat.p,
                                 h_mat.row,
                                 psi_in.data<T>(),
                                 1,
                                 zero_,
                                 hpsi_out.data<T>(),
                                 1);

            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
        };

        auto spsi_func = [s_mat, one_, zero_](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");

            Device* ctx = {};

            gemv_op<T, Device>()(ctx,
                                 'N',
                                 s_mat.row,
                                 s_mat.col,
                                 one_,
                                 s_mat.p,
                                 s_mat.row,
                                 psi_in.data<T>(),
                                 1,
                                 zero_,
                                 spsi_out.data<T>(),
                                 1);

            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
        };

        // if (this->is_first_scf)
        // {
        for (size_t i = 0; i < psi.get_nbands(); i++)
        {
            for (size_t j = 0; j < psi.get_nbasis(); j++)
            {
                psi(i, j) = *zero_;
            }
            psi(i, i) = *one_;
        }
        // }

        auto psi_tensor = ct::TensorMap(psi.get_pointer(),
                                        ct::DataTypeToEnum<T>::value,
                                        ct::DeviceTypeToEnum<ct_Device>::value,
                                        ct::TensorShape({psi.get_nbands(), psi.get_nbasis()}))
                              .slice({0, 0}, {psi.get_nbands(), psi.get_current_nbas()});

        auto eigen_tensor = ct::TensorMap(eigenvalue,
                                          ct::DataTypeToEnum<Real>::value,
                                          ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                          ct::TensorShape({psi.get_nbands()}));

        auto prec_tensor = ct::TensorMap(this->precondition_lcao.data(),
                                         ct::DataTypeToEnum<Real>::value,
                                         ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                         ct::TensorShape({static_cast<int>(this->precondition_lcao.size())}))
                               .slice({0}, {psi.get_current_nbas()});

        cg.diag(hpsi_func, spsi_func, psi_tensor, eigen_tensor, prec_tensor);

        // TODO: Double check tensormap's potential problem
        ct::TensorMap(psi.get_pointer(), psi_tensor, {psi.get_nbands(), psi.get_nbasis()}).sync(psi_tensor);
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
        auto psi_pool = psi::Psi<T>(1, ncol_bands_pool, k2d.get_p2D_pool()->nrow, nullptr);
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

template class HSolverLCAO<double>;
template class HSolverLCAO<std::complex<double>>;

} // namespace hsolver
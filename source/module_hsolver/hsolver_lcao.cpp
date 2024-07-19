#include "hsolver_lcao.h"

#include "diago_cg.h"

#ifdef __MPI
#include "diago_scalapack.h"
#else
#include "diago_lapack.h"
#endif

#include "module_base/timer.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/write_HS.h"

#include "module_base/global_variable.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_types.h>
#ifdef __CUSOLVERMP
#include "diago_cusolvermp.h"
#endif // __CUSOLVERMP
#ifdef __ELPA
#include "diago_elpa.h"
#endif
#ifdef __CUDA
#include "diago_cusolver.h"
#endif
#ifdef __PEXSI
#include "diago_pexsi.h"
#include "module_elecstate/elecstate_lcao.h"
#endif

namespace hsolver
{

template <typename T, typename Device>
void HSolverLCAO<T, Device>::solve(hamilt::Hamilt<T>* pHamilt,
                                   psi::Psi<T>& psi,
                                   elecstate::ElecState* pes,
                                   const std::string method_in,
                                   const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLCAO", "solve");
    ModuleBase::timer::tick("HSolverLCAO", "solve");
    // select the method of diagonalization
    this->method = method_in;

#ifdef __PEXSI
    if (this->method == "pexsi")
    {
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete[] this->pdiagh;
                this->pdiagh = nullptr;
            }
            auto tem = dynamic_cast<DiagoPexsi<T>*>(this->pdiagh);
        }
        if (this->pdiagh == nullptr)
        {
            DiagoPexsi<T>* tem = new DiagoPexsi<T>(this->ParaV);
            this->pdiagh = tem;
            // this->pdiagh = dynamic_cast<DiagoPexsi<T>*>(tem);
            this->pdiagh->method = this->method;
        }
    }
#endif


    // Zhang Xiaoyang :  Please modify Pesxi usage later
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

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);

        // solve eigenvector and eigenvalue for H(k)
        this->hamiltSolvePsiK(pHamilt, psi, &(pes->ekb(ik, 0)));
    }

    if (this->method == "cg_in_lcao")
    {
        this->is_first_scf = false;
    }

    if (this->method != "genelpa" && this->method != "scalapack_gvx" && this->method != "lapack"
        && this->method != "cusolver" && this->method != "cusolvermp" && this->method != "cg_in_lcao"
        && this->method != "pexsi")
    {
        //delete this->pdiagh;
        //this->pdiagh = nullptr;
    }

    // used in nscf calculation
    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverLCAO", "solve");
        return;
    }

    // calculate charge by psi
    // called in scf calculation
#ifdef __PEXSI
    if (this->method == "pexsi")
    {
        DiagoPexsi<T> tem = dynamic_cast<DiagoPexsi<T>*>(this->pdiagh);
        if (tem == nullptr)
            ModuleBase::WARNING_QUIT("HSolverLCAO", "pexsi need debug!");
        elecstate::ElecStateLCAO<T>* _pes = dynamic_cast<elecstate::ElecStateLCAO<T>*>(pes);
        pes->f_en.eband = tem->totalFreeEnergy;
        // maybe eferm could be dealt with in the future
        _pes->dmToRho(tem->DM, tem->EDM);
    }
    else
#endif
    {
        pes->psiToRho(psi);
    }
    ModuleBase::timer::tick("HSolverLCAO", "solve");
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
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");
#endif
#ifdef __PEXSI
    else if (this->method == "pexsi")
    {
        this->pdiagh->diag(hm, psi, eigenvalue);
    }
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

        DiagoCG<T, Device> cg(GlobalV::BASIS_TYPE,
                              GlobalV::CALCULATION,
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
                h_mat.p[h_mat.row * j + i] = hsolver::my_conj(h_mat.p[h_mat.row * i + j]);
                s_mat.p[s_mat.row * j + i] = hsolver::my_conj(s_mat.p[s_mat.row * i + j]);
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

        if (this->is_first_scf)
        {
            for (size_t i = 0; i < psi.get_nbands(); i++)
            {
                for (size_t j = 0; j < psi.get_nbasis(); j++)
                {
                    psi(i, j) = *zero_;
                }
                psi(i, i) = *one_;
            }
        }

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

template class HSolverLCAO<double>;
template class HSolverLCAO<std::complex<double>>;

} // namespace hsolver
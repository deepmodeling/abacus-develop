#include "diago_iter_assist.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/device.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/kernels/hegvd_op.h"

#include <cassert>

namespace hsolver
{

//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart vectors psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace(const HSOperator<T, Device>& op,
                                               const T* psi,
                                               T* evc,
                                               const int nstart,
                                               const int n_band,
                                               const int dmin,
                                               const int dmax,
                                               Real* en,
                                               const diag_comm_info& diag_comm,
                                               const bool S_orth)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_subspace");
    ModuleBase::timer::start("DiagoIterAssist", "diag_subspace");

    assert(n_band <= nstart);

    // scc is overlap (optional, only needed if input is not s-orthogonal)
    T *hcc = nullptr, *scc = nullptr, *vcc = nullptr;

    // hcc is reduced hamiltonian matrix
    resmem_complex_op()(hcc, nstart * nstart, "DiagSub::hcc");
    setmem_complex_op()(hcc, 0, nstart * nstart);

    // scc is overlap matrix, only needed when psi is not orthogonal
    if(!S_orth){
        resmem_complex_op()(scc, nstart * nstart, "DiagSub::scc");
        setmem_complex_op()(scc, 0, nstart * nstart);
    }
    
    // vcc is eigenvector matrix of the reduced generalized eigenvalue problem
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(vcc, 0, nstart * nstart);

    // temp holds H|psi>, then S|psi>, then the rotated vectors; it is separate
    // from evc so that evc may alias psi
    T* temp = nullptr;
    resmem_complex_op()(temp, nstart * dmax, "DiagSub::temp");
    setmem_complex_op()(temp, 0, nstart * dmax);

    { // code block to calculate hcc and scc
        T *hpsi = temp;
        op.hpsi(psi, hpsi, dmax, nstart);

        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         nstart,
                                         nstart,
                                         dmin,
                                         &one,
                                         psi,
                                         dmax,
                                         hpsi,
                                         dmax,
                                         &zero,
                                         hcc,
                                         nstart);

        if(!S_orth){
            // Only calculate S_sub if not orthogonal
            T *spsi = temp;
            op.spsi(psi, spsi, dmax, nstart);

            ModuleBase::gemm_op<T, Device>()('C',
                                            'N',
                                            nstart,
                                            nstart,
                                            dmin,
                                            &one,
                                            psi,
                                            dmax,
                                            spsi,
                                            dmax,
                                            &zero,
                                            scc,
                                            nstart);
        }
    }

    if (diag_comm.nproc > 1)
    {
#ifdef __MPI
        Parallel_Common::reduce_dev<T, Device>(hcc, nstart * nstart, diag_comm.comm);
        if(!S_orth){
            Parallel_Common::reduce_dev<T, Device>(scc, nstart * nstart, diag_comm.comm);
        }
#endif
    }

    // after generation of H and (optionally) S matrix, diag them
    if (S_orth) {
        // Solve standard eigenproblem: H_sub * y = lambda * y
        DiagoIterAssist::diag_heevx(nstart, n_band, hcc, nstart, en, vcc);
    } else {
        // Solve generalized eigenproblem: H_sub * y = lambda * S_sub * y
        DiagoIterAssist::diag_hegvd(nstart, n_band, hcc, scc, nstart, en, vcc);
    }

    { // code block to calculate evc
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         dmin,
                                         n_band,
                                         nstart,
                                         &one,
                                         psi, // dmin * nstart
                                         dmax,
                                         vcc, // nstart * n_band
                                         nstart,
                                         &zero,
                                         temp,
                                         dmin);
    }

    ModuleBase::matrixCopy<T, Device>()(n_band, dmin, temp, dmin, evc, dmax);

    delmem_complex_op()(temp);
    delmem_complex_op()(hcc);
    if(!S_orth){
        delmem_complex_op()(scc);
    }
    delmem_complex_op()(vcc);

    ModuleBase::timer::end("DiagoIterAssist", "diag_subspace");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace(const HSOperator<T, Device>& op,
                                               const psi::Psi<T, Device>& psi, // [in] wavefunction
                                               psi::Psi<T, Device>& evc,       // [out] wavefunction, eigenvectors
                                               Real* en,                       // [out] eigenvalues
                                               const diag_comm_info& diag_comm,
                                               int n_band,       // [in] number of bands to be calculated, also number of rows
                                                                 // of evc, if set to 0, n_band = nstart, default 0
                                               const bool S_orth // [in] if true, psi is assumed to be already S-orthogonalized
)
{
    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)
    const int nstart = psi.get_nbands();
    // n_band = 0 means default, set n_band = nstart
    if (n_band == 0)
    {
        n_band = nstart;
    }

    // dmin is the active number of plane waves or atomic orbitals
    // dmax is the leading dimension of psi
    diag_subspace(op,
                  psi.get_pointer(),
                  evc.get_pointer(),
                  nstart,
                  n_band,
                  psi.get_current_ngk(),
                  psi.get_nbasis(),
                  en,
                  diag_comm,
                  S_orth);
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace_init(const HSOperator<T, Device>& op,
                                                    const T* psi,
                                                    int psi_nr,
                                                    int psi_nc,
                                                    psi::Psi<T, Device>& evc,
                                                    Real* en,
                                                    const std::string& basis_type,
                                                    const std::string& calculation,
                                                    const diag_comm_info& diag_comm)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_subspace_init");
    ModuleBase::timer::start("DiagoIterAssist", "diag_subspace_init");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)

    const int nstart = psi_nr;
    const int n_band = evc.get_nbands();
    const int dmax = evc.get_nbasis();
    const int dmin = evc.get_current_ngk();

    T *hcc = nullptr, *scc = nullptr, *vcc = nullptr;
    resmem_complex_op()(hcc, nstart * nstart, "DiagSub::hcc");
    resmem_complex_op()(scc, nstart * nstart, "DiagSub::scc");
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(hcc, 0, nstart * nstart);
    setmem_complex_op()(scc, 0, nstart * nstart);
    setmem_complex_op()(vcc, 0, nstart * nstart);

    if (base_device::get_device_type(ctx) == base_device::GpuDevice)
    {
        // band by band on the GPU: the scratch buffer holds one vector only
        T* temp = nullptr;
        resmem_complex_op()(temp, psi_nc, "DiagSub::temp");
        setmem_complex_op()(temp, 0, psi_nc);

        T* hpsi = temp;
        for (int i = 0; i < nstart; i++)
        {
            // H|Psi> to get hpsi for target band
            op.hpsi(psi + i * psi_nc, hpsi, psi_nc, 1);

            // calculate the related elements in hcc <Psi|H|Psi>
            ModuleBase::gemv_op<T, Device>()('C', psi_nc, nstart, &one, psi, psi_nc, hpsi, 1, &zero, hcc + i * nstart, 1);
        }

        T* spsi = temp;
        for (int i = 0; i < nstart; i++)
        {
            op.spsi(psi + i * psi_nc, spsi, psi_nc, 1);

            ModuleBase::gemv_op<T, Device>()('C',
                                             psi_nc,
                                             nstart,
                                             &one,
                                             psi,
                                             psi_nc, // nbasis
                                             spsi,
                                             1,
                                             &zero,
                                             scc + i * nstart,
                                             1);
        }
        delmem_complex_op()(temp);
    }
    else if (base_device::get_device_type(ctx) == base_device::CpuDevice)
    {
        // hpsi and spsi share the temp space
        T* temp = nullptr;
        resmem_complex_op()(temp, nstart * psi_nc, "DiagSub::temp");
        setmem_complex_op()(temp, 0, nstart * psi_nc);

        T* hpsi = temp;
        op.hpsi(psi, hpsi, psi_nc, nstart);

        ModuleBase::gemm_op<T, Device>()('C', 'N', nstart, nstart, dmin, &one, psi, psi_nc, hpsi, psi_nc, &zero, hcc, nstart);

        T* spsi = temp;
        op.spsi(psi, spsi, psi_nc, nstart);

        ModuleBase::gemm_op<T, Device>()('C', 'N', nstart, nstart, dmin, &one, psi, psi_nc, spsi, psi_nc, &zero, scc, nstart);
        delmem_complex_op()(temp);
    }

    // a Hamiltonian may carry a term hpsi() does not cover (EXX in lcao_in_pw)
    op.add_to_subspace_h(hcc, nstart);

    if (diag_comm.nproc > 1)
    {
#ifdef __MPI
        Parallel_Common::reduce_dev<T, Device>(hcc, nstart * nstart, diag_comm.comm);
        Parallel_Common::reduce_dev<T, Device>(scc, nstart * nstart, diag_comm.comm);
#endif
    }

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diag_hegvd(nstart, n_band, hcc, scc, nstart, en, vcc);

    op.export_subspace_vec(vcc, nstart, n_band);

    //=======================
    // diagonize the H-matrix
    //=======================
    if ((basis_type == "lcao" || basis_type == "lcao_in_pw") && calculation == "nscf")
    {
        // The caller requested eigenvalues only, so no wavefunction rotation is needed.
    }
    else if ((basis_type == "lcao" || basis_type == "lcao_in_pw" || basis_type == "pw")
             && (calculation == "scf" || calculation == "md"
                 || calculation == "relax")) // pengfei 2014-10-13
    {
        // because psi and evc are different here,
        // I think if psi and evc are the same,
        // there may be problems, mohan 2011-01-01
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         dmax,
                                         n_band,
                                         nstart,
                                         &one,
                                         psi, // dmax * nstart
                                         dmax,
                                         vcc, // nstart * n_band
                                         nstart,
                                         &zero,
                                         evc.get_pointer(),
                                         dmax);
    }
    else
    {
        assert(psi != evc.get_pointer());

        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         dmin,
                                         n_band,
                                         nstart,
                                         &one,
                                         psi, // dmin * nstart
                                         dmax,
                                         vcc, // nstart * n_band
                                         nstart,
                                         &zero,
                                         evc.get_pointer(),
                                         dmax);
    }

    delmem_complex_op()(hcc);
    delmem_complex_op()(scc);
    delmem_complex_op()(vcc);
    ModuleBase::timer::end("DiagoIterAssist", "diag_subspace_init");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_heevx(const int matrix_size,
                                                       const int num_eigenpairs,
                                                       const T *h,
                                                       const int ldh,
                                                       Real *e, // always in CPU
                                                       T *v)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_heevx");
    ModuleBase::timer::start("DiagoIterAssist", "diag_heevx");

    Real *eigenvalues = nullptr;
    // device memory for eigenvalues
    resmem_var_op()(eigenvalues, matrix_size);
    setmem_var_op()(eigenvalues, 0, matrix_size);

    // (const Device *d, const int matrix_size, const int lda, const T *A, const int num_eigenpairs, Real *eigenvalues, T *eigenvectors);
    heevx_op<T, Device>()(ctx, matrix_size, ldh, h, num_eigenpairs, eigenvalues, v);

    if (base_device::get_device_type(ctx) == base_device::GpuDevice)
    {
#if ((defined __CUDA) || (defined __ROCM))
        // eigenvalues to e, from device to host
        syncmem_var_d2h_op()(e, eigenvalues, num_eigenpairs);
#endif
    }
    else if (base_device::get_device_type(ctx) == base_device::CpuDevice)
    {
        // eigenvalues to e
        syncmem_var_op()(e, eigenvalues, num_eigenpairs);
    }

    delmem_var_op()(eigenvalues);

    ModuleBase::timer::end("DiagoIterAssist", "diag_heevx");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_hegvd(const int nstart,
                                              const int nbands,
                                              const T *hcc,
                                              T *scc,
                                              const int ldh, // nstart
                                              Real *e,       // always in CPU
                                              T *vcc)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_hegvd");
    ModuleBase::timer::start("DiagoIterAssist", "diag_hegvd");

    Real *eigenvalues = nullptr;
    resmem_var_op()(eigenvalues, nstart);
    setmem_var_op()(eigenvalues, 0, nstart);

    hegvd_op<T, Device>()(ctx, nstart, ldh, hcc, scc, eigenvalues, vcc);

    if (base_device::get_device_type(ctx) == base_device::GpuDevice)
    {
#if ((defined __CUDA) || (defined __ROCM))
        // set eigenvalues in GPU to e in CPU
        syncmem_var_d2h_op()(e, eigenvalues, nbands);
#endif
    }
    else if (base_device::get_device_type(ctx) == base_device::CpuDevice)
    {
        // set eigenvalues in CPU to e in CPU
        syncmem_var_op()(e, eigenvalues, nbands);
    }

    delmem_var_op()(eigenvalues);

    // const bool all_eigenvalues = (nstart == nbands);
    // if (all_eigenvalues) {
    //     //===========================
    //     // calculate all eigenvalues
    //     //===========================
    //     // dngv_op<Real, Device>()(ctx, nstart, ldh, hcc, scc, res, vcc);
    //     dngvd_op<Real, Device>()(ctx, nstart, ldh, hcc, scc, res, vcc);
    // }
    // else {
    //     //=====================================
    //     // calculate only m lowest eigenvalues
    //     //=====================================
    //     dngvx_op<Real, Device>()(ctx, nstart, ldh, hcc, scc, nbands, res, vcc);
    // }

    ModuleBase::timer::end("DiagoIterAssist", "diag_hegvd");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::cal_hs_subspace(const HSOperator<T, Device>& op,
                                                 const psi::Psi<T, Device>& psi, // [in] wavefunction
                                                 T* hcc,
                                                 T* scc,
                                                 const diag_comm_info& diag_comm)
{
    const int nstart = psi.get_nbands();
    
    setmem_complex_op()(hcc, 0, nstart * nstart);
    setmem_complex_op()(scc, 0, nstart * nstart);

    const int dmin = psi.get_current_ngk();
    const int dmax = psi.get_nbasis();

    T* temp = nullptr;
    resmem_complex_op()(temp, nstart * dmax, "DiagSub::temp");
    setmem_complex_op()(temp, 0, nstart * dmax);

    { // code block to calculate hcc and scc
        T* hpsi = temp;
        op.hpsi(psi.get_pointer(), hpsi, dmax, nstart);

        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         nstart,
                                         nstart,
                                         dmin,
                                         &one,
                                         psi.get_pointer(),
                                         dmax,
                                         hpsi,
                                         dmax,
                                         &zero,
                                         hcc,
                                         nstart);

        T* spsi = temp;
        op.spsi(psi.get_pointer(), spsi, dmax, nstart);

        ModuleBase::gemm_op<T, Device>()('C',
                                         'N',
                                         nstart,
                                         nstart,
                                         dmin,
                                         &one,
                                         psi.get_pointer(),
                                         dmax,
                                         spsi,
                                         dmax,
                                         &zero,
                                         scc,
                                         nstart);
    }

    if (diag_comm.nproc > 1)
    {
#ifdef __MPI
        Parallel_Common::reduce_dev<T, Device>(hcc, nstart * nstart, diag_comm.comm);
        Parallel_Common::reduce_dev<T, Device>(scc, nstart * nstart, diag_comm.comm);
#endif
    }

    delmem_complex_op()(temp);
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_responce( const T* hcc,
                                                T* scc,
                                                const int nbands,
                                                const T* mat_in,           // [out] target matrix to be multiplied
                                                T* mat_out,
                                                int mat_col,          // [in] number of columns of target matrix
                                                Real* en                           // [out] eigenvalues
)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_responce");
    ModuleBase::timer::start("DiagoIterAssist", "diag_responce");

    const int nstart = nbands;

    T *vcc = nullptr;
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(vcc, 0, nstart * nstart);

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diag_hegvd(nstart, nstart, hcc, scc, nstart, en, vcc);

    { // code block to calculate tar_mat
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         mat_col,
                                         nstart,
                                         nstart,
                                         &one,
                                         mat_in, // mat_col * nstart
                                         mat_col,
                                         vcc, // nstart * nstart
                                         nstart,
                                         &zero,
                                         mat_out,
                                         mat_col);
    }

    delmem_complex_op()(vcc);

    ModuleBase::timer::end("DiagoIterAssist", "diag_responce");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace_psi(const T* hcc,
                              T* scc,
                              const int dim_subspace,
                              psi::Psi<T, Device>& evc,
                              Real* en
)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_subspace_psi");
    ModuleBase::timer::start("DiagoIterAssist", "diag_subspace_psi");

    const int nstart = dim_subspace;
    const int n_band = evc.get_nbands();

    T *vcc = nullptr;
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(vcc, 0, nstart * nstart);

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diag_hegvd(nstart, nstart, hcc, scc, nstart, en, vcc);

    { // code block to calculate tar_mat
        const int dmin = evc.get_current_ngk();
        const int dmax = evc.get_nbasis();
        T* temp = nullptr;
        resmem_complex_op()(temp, nstart * dmax, "DiagSub::temp");
        setmem_complex_op()(temp, 0, nstart * dmax);
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         dmin,
                                         n_band,
                                         nstart,
                                         &one,
                                         evc.get_pointer(), // dmin * nstart
                                         dmax,
                                         vcc, // nstart * n_band
                                         nstart,
                                         &zero,
                                         temp,
                                         dmin);
        ModuleBase::matrixCopy<T, Device>()(n_band, dmin, temp, dmin, evc.get_pointer(), dmax);
        delmem_complex_op()(temp);
    }

    delmem_complex_op()(vcc);

    ModuleBase::timer::end("DiagoIterAssist", "diag_subspace_psi");
}

template class DiagoIterAssist<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoIterAssist<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>;
#endif

#ifdef __LCAO
template class DiagoIterAssist<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoIterAssist<double, base_device::DEVICE_GPU>;
#endif
#endif
} // namespace hsolver

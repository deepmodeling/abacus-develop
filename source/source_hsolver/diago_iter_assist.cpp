#include "diago_iter_assist.h"
#include "source_io/module_parameter/parameter.h"
#include "source_base/complexmatrix.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/module_device/device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_base/kernels/math_kernel_op.h"

namespace hsolver
{

//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace(const hamilt::Hamilt<T, Device>* const pHamilt, // hamiltonian operator carrier
                                                const psi::Psi<T, Device>& psi,     // [in] wavefunction
                                                psi::Psi<T, Device>& evc,           // [out] wavefunction, eigenvectors
                                                Real* en,                           // [out] eigenvalues
                                                int n_band, // [in] number of bands to be calculated, also number of rows
                                                           // of evc, if set to 0, n_band = nstart, default 0
                                                const bool S_orth // [in] if true, psi is assumed to be already S-orthogonalized
)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_subspace");
    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)
    const int nstart = psi.get_nbands();
    // n_band = 0 means default, set n_band = nstart
    if (n_band == 0)
    {
        n_band = nstart;
    }
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

    // dmin is the active number of plane waves or atomic orbitals
    // dmax is the leading dimension of psi
    const int dmin = psi.get_current_ngk();
    const int dmax = psi.get_nbasis();

    T *temp = nullptr; /// temporary array for calculation of evc
    bool in_place = false; ///< if temp and evc share the same memory
    if (psi.get_pointer() != evc.get_pointer() && psi.get_nbands() == evc.get_nbands())
    { // use memory of evc as temp
        temp = evc.get_pointer();
        in_place = true;
    }
    else
    {
        resmem_complex_op()(temp, nstart * dmax, "DiagSub::temp");
    }

    { // code block to calculate hcc and scc
        setmem_complex_op()(temp, 0, nstart * dmax);

        T *hpsi = temp;
        // do hPsi for all bands
        psi::Range all_bands_range(1, psi.get_current_k(), 0, nstart - 1);
        hpsi_info hpsi_in(&psi, all_bands_range, hpsi);
        pHamilt->ops->hPsi(hpsi_in);

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

        if(!S_orth){
            // Only calculate S_sub if not orthogonal
            T *spsi = temp;
            // do sPsi for all bands
            pHamilt->sPsi(psi.get_pointer(), spsi, dmax, dmin, nstart);

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
    }

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_pool(hcc, nstart * nstart);
        if(!S_orth){
            Parallel_Reduce::reduce_pool(scc, nstart * nstart);
        }
    }

    // after generation of H and (optionally) S matrix, diag them
    if (S_orth) {
        // Solve standard eigenproblem: H_sub * y = lambda * y
        DiagoIterAssist::diag_heevx(nstart, n_band, hcc, nstart, en, vcc);
    } else {
        // Solve generalized eigenproblem: H_sub * y = lambda * S_sub * y
        DiagoIterAssist::diag_hegvd(nstart, n_band, hcc, scc, nstart, en, vcc);
    }

    const int ld_temp = in_place ? dmax : dmin;

    { // code block to calculate evc
        ModuleBase::gemm_op<T, Device>()('N',
                                         'N',
                                         dmin,
                                         n_band,
                                         nstart,
                                         &one,
                                         psi.get_pointer(), // dmin * nstart
                                         dmax,
                                         vcc, // nstart * n_band
                                         nstart,
                                         &zero,
                                         temp,
                                         ld_temp);
    }

    if (!in_place)
    {
        ModuleBase::matrixCopy<T, Device>()(n_band, ld_temp, temp, ld_temp, evc.get_pointer(), dmax);
        delmem_complex_op()(temp);
    }
    delmem_complex_op()(hcc);
    if(!S_orth){
        delmem_complex_op()(scc);
    }
    delmem_complex_op()(vcc);

    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diag_subspace_init(hamilt::Hamilt<T, Device>* pHamilt,
    const T* psi,
    int psi_nr,
    int psi_nc,
    psi::Psi<T, Device>& evc,
    Real* en,
    const std::function<void(T*, const int)>& add_to_hcc,
    const std::function<void(const T* const, const int, const int)>& export_vcc)
{
    ModuleBase::TITLE("DiagoIterAssist", "diag_subspace_init");
    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace_init");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)

    const int nstart = psi_nr;
    const int n_band = evc.get_nbands();
    const int dmax = evc.get_nbasis();
    const int dmin = evc.get_current_ngk();

    // skip the diagonalization if the operators are not allocated
    if (pHamilt->ops == nullptr)
    {
        ModuleBase::WARNING(
            "DiagoIterAssist::diag_subspace_init",
            "Severe warning: Operators in Hamilt are not allocated yet, will return value of psi to evc directly\n");
        for (int iband = 0; iband < n_band; iband++)
        {
            for (int ig = 0; ig < dmax; ig++)
            {
                evc(iband, ig) = psi[iband * dmax + ig];
            }
            en[iband] = 0.0;
        }
        return;
    }

    // ModuleBase::ComplexMatrix hc(nstart, nstart);
    // ModuleBase::ComplexMatrix sc(nstart, nstart);
    // ModuleBase::ComplexMatrix hvec(nstart, n_band);
    T *hcc = nullptr, *scc = nullptr, *vcc = nullptr;
    resmem_complex_op()(hcc, nstart * nstart, "DiagSub::hcc");
    resmem_complex_op()(scc, nstart * nstart, "DiagSub::scc");
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(hcc, 0, nstart * nstart);
    setmem_complex_op()(scc, 0, nstart * nstart);
    setmem_complex_op()(vcc, 0, nstart * nstart);

    if (base_device::get_device_type(ctx) == base_device::GpuDevice)
    {
        psi::Psi<T, Device> psi_temp(1, 1, psi_nc, dmin, true);

        T* ppsi = psi_temp.get_pointer();
        // hpsi and spsi share the temp space
        T* temp = nullptr;
        resmem_complex_op()(temp, psi_nc, "DiagSub::temp");
        setmem_complex_op()(temp, 0, psi_nc);

        T* hpsi = temp;
        // do hPsi band by band
        for (int i = 0; i < nstart; i++)
        {
            // psi_temp is one band psi, psi is all bands psi, the range always is 1 for the only band in psi_temp
            syncmem_complex_op()(ppsi, psi + i * psi_nc, psi_nc);
            psi::Range band_by_band_range(true, 0, 0, 0);
            hpsi_info hpsi_in(&psi_temp, band_by_band_range, hpsi);

            // H|Psi> to get hpsi for target band
            pHamilt->ops->hPsi(hpsi_in);

            // calculate the related elements in hcc <Psi|H|Psi>
            ModuleBase::gemv_op<T, Device>()('C', psi_nc, nstart, &one, psi, psi_nc, hpsi, 1, &zero, hcc + i * nstart, 1);
        }

        T* spsi = temp;
        // do sPsi band by band
        for (int i = 0; i < nstart; i++)
        {
            syncmem_complex_op()(ppsi, psi + i * psi_nc, psi_nc);
            pHamilt->sPsi(ppsi, spsi, dmin, dmin, 1);

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
        psi::Psi<T, Device> psi_temp(1, nstart, psi_nc, dmin, true);

        T* ppsi = psi_temp.get_pointer();
        syncmem_complex_op()(ppsi, psi, psi_temp.size());
        // hpsi and spsi share the temp space
        T* temp = nullptr;
        resmem_complex_op()(temp, nstart * psi_nc, "DiagSub::temp");
        setmem_complex_op()(temp, 0, nstart * psi_nc);

        T* hpsi = temp;
        // do hPsi for all bands
        psi::Range all_bands_range(true, 0, 0, nstart - 1);
        hpsi_info hpsi_in(&psi_temp, all_bands_range, hpsi);
        pHamilt->ops->hPsi(hpsi_in);

        ModuleBase::gemm_op<T, Device>()('C', 'N', nstart, nstart, dmin, &one, ppsi, dmax, hpsi, dmax, &zero, hcc, nstart);

        T* spsi = temp;
        // do sPsi for all bands
        pHamilt->sPsi(ppsi, spsi, psi_temp.get_nbasis(), psi_temp.get_nbasis(), psi_temp.get_nbands());

        ModuleBase::gemm_op<T, Device>()('C', 'N', nstart, nstart, dmin, &one, ppsi, dmax, spsi, dmax, &zero, scc, nstart);
        delmem_complex_op()(temp);

        add_to_hcc(hcc, nstart);

    }

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_pool(hcc, nstart * nstart);
        Parallel_Reduce::reduce_pool(scc, nstart * nstart);
    }

    // after generation of H and S matrix, diag them
    /// this part only for test, eigenvector would have different phase caused by micro numerical perturbation
    /// set 8 bit effective accuracy would help for debugging
    /*for(int i=0;i<nstart;i++)
    {
        for(int j=0;j<nstart;j++)
        {
            if(std::norm(hc(i,j))<1e-10) hc(i,j) = ModuleBase::ZERO;
            else hc(i,j) = std::complex<double>(double(int(hc(i,j).real()*100000000))/100000000, 0);
            if(std::norm(sc(i,j))<1e-10) sc(i,j) = ModuleBase::ZERO;
            else sc(i,j) = std::complex<double>(double(int(sc(i,j).real()*100000000))/100000000, 0);
        }
    }*/

    DiagoIterAssist::diag_hegvd(nstart, n_band, hcc, scc, nstart, en, vcc);

    export_vcc(vcc, nstart, n_band);

    //=======================
    // diagonize the H-matrix
    //=======================
    if ((PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw") && PARAM.inp.calculation == "nscf")
    {
        GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
    }
    else if ((PARAM.inp.basis_type == "lcao" || PARAM.inp.basis_type == "lcao_in_pw" || PARAM.inp.basis_type == "pw")
             && (PARAM.inp.calculation == "scf" || PARAM.inp.calculation == "md"
                 || PARAM.inp.calculation == "relax")) // pengfei 2014-10-13
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

        // T* evctemp = nullptr;
        // resmem_complex_op()(ctx, evctemp, n_band * dmin, "DiagSub::evctemp");
        // setmem_complex_op()(ctx, evctemp, 0, n_band * dmin);

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

        // matrixCopy<T, Device>()(ctx, n_band, evctemp, dmin, evc.get_pointer(), dmax);

        // delmem_complex_op()(ctx, evctemp);
    }

    delmem_complex_op()(hcc);
    delmem_complex_op()(scc);
    delmem_complex_op()(vcc);
    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace_init");
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
    ModuleBase::timer::tick("DiagoIterAssist", "diag_heevx");

    Real *eigenvalues = nullptr;
    // device memory for eigenvalues
    resmem_var_op()(eigenvalues, matrix_size);
    setmem_var_op()(eigenvalues, 0, matrix_size);

    // (const Device *d, const int matrix_size, const int lda, const T *A, const int num_eigenpairs, Real *eigenvalues, T *eigenvectors);
    heevx_op<T, Device>()(ctx, matrix_size, ldh, h, num_eigenpairs, eigenvalues, v);

    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
#if ((defined __CUDA) || (defined __ROCM))
        // eigenvalues to e, from device to host
        syncmem_var_d2h_op()(e, eigenvalues, num_eigenpairs);
#endif
    }
    else if (base_device::get_device_type<Device>(ctx) == base_device::CpuDevice)
    {
        // eigenvalues to e
        syncmem_var_op()(e, eigenvalues, num_eigenpairs);
    }

    delmem_var_op()(eigenvalues);

    ModuleBase::timer::tick("DiagoIterAssist", "diag_heevx");
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
    ModuleBase::timer::tick("DiagoIterAssist", "diag_hegvd");

    Real *eigenvalues = nullptr;
    resmem_var_op()(eigenvalues, nstart);
    setmem_var_op()(eigenvalues, 0, nstart);

    hegvd_op<T, Device>()(ctx, nstart, ldh, hcc, scc, eigenvalues, vcc);

    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
#if ((defined __CUDA) || (defined __ROCM))
        // set eigenvalues in GPU to e in CPU
        syncmem_var_d2h_op()(e, eigenvalues, nbands);
#endif
    }
    else if (base_device::get_device_type<Device>(ctx) == base_device::CpuDevice)
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

    ModuleBase::timer::tick("DiagoIterAssist", "diag_hegvd");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::cal_hs_subspace(const hamilt::Hamilt<T, Device>* pHamilt, // hamiltonian operator carrier
                                                const psi::Psi<T, Device>& psi,     // [in] wavefunction
                                                T *hcc, 
                                                T *scc)
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
        setmem_complex_op()(temp, 0, nstart * dmax);

        T* hpsi = temp;
        // do hPsi for all bands
        psi::Range all_bands_range(1, psi.get_current_k(), 0, nstart - 1);
        hpsi_info hpsi_in(&psi, all_bands_range, hpsi);
        pHamilt->ops->hPsi(hpsi_in);

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
        // do sPsi for all bands
        pHamilt->sPsi(psi.get_pointer(), spsi, dmax, dmin, nstart);

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

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_pool(hcc, nstart * nstart);
        Parallel_Reduce::reduce_pool(scc, nstart * nstart);
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
    ModuleBase::timer::tick("DiagoIterAssist", "diag_responce");

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

    ModuleBase::timer::tick("DiagoIterAssist", "diag_responce");
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
    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace_psi");

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

    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace_psi");
}

template <typename T, typename Device>
bool DiagoIterAssist<T, Device>::test_exit_cond(const int& ntry, const int& notconv)
{
    //================================================================
    // If this logical function is true, need to do diag_subspace
    // and cg again.
    //================================================================

    bool scf = true;
    if (PARAM.inp.calculation == "nscf") {
        scf = false;
}

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ((!scf && (notconv > 0)));

    // if self consistent calculation, if not converged > 5,
    // using diag_subspace and cg method again. ntry++
    const bool f3 = ((scf && (notconv > 5)));
    return (f1 && (f2 || f3));
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

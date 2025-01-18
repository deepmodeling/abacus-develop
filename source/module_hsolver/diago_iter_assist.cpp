#include "diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "module_base/blas_connector.h"
#include "module_base/complexmatrix.h"
#include "module_base/constants.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/module_device/device.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

namespace hsolver
{

//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diagH_subspace(const hamilt::Hamilt<T, Device>* const pHamilt, // hamiltonian operator carrier
                                                const psi::Psi<T, Device>& psi,     // [in] wavefunction
                                                psi::Psi<T, Device>& evc,           // [out] wavefunction
                                                Real* en,                           // [out] eigenvalues
                                                int n_band // [in] number of bands to be calculated, also number of rows
                                                           // of evc, if set to 0, n_band = nstart, default 0
)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_subspace");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)
    const int nstart = psi.get_nbands();
    if (n_band == 0)
    {
        n_band = nstart;
    }
    assert(n_band <= nstart);

    T *hcc = nullptr, *scc = nullptr, *vcc = nullptr;
    resmem_complex_op()(hcc, nstart * nstart, "DiagSub::hcc");
    resmem_complex_op()(scc, nstart * nstart, "DiagSub::scc");
    resmem_complex_op()(vcc, nstart * nstart, "DiagSub::vcc");
    setmem_complex_op()(hcc, 0, nstart * nstart);
    setmem_complex_op()(scc, 0, nstart * nstart);
    setmem_complex_op()(vcc, 0, nstart * nstart);

    const int dmin = psi.get_current_ngk();
    const int dmax = psi.get_nbasis();

    T* temp = nullptr;
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

        T* hphi = temp;
        // do hPsi for all bands
        psi::Range all_bands_range(1, psi.get_current_k(), 0, nstart - 1);
        hpsi_info hpsi_in(&psi, all_bands_range, hphi);
        pHamilt->ops->hPsi(hpsi_in);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             &one,
                             psi.get_pointer(),
                             dmax,
                             hphi,
                             dmax,
                             &zero,
                             hcc,
                             nstart);

        T* sphi = temp;
        // do sPsi for all bands
        pHamilt->sPsi(psi.get_pointer(), sphi, dmax, dmin, nstart);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             &one,
                             psi.get_pointer(),
                             dmax,
                             sphi,
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

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diagH_LAPACK(nstart, n_band, hcc, scc, nstart, en, vcc);


    const int ld_temp = in_place ? dmax : dmin;

    { // code block to calculate evc
        gemm_op<T, Device>()(ctx,
                             'N',
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
        matrixSetToAnother<T, Device>()(ctx, n_band, temp, ld_temp, evc.get_pointer(), dmax);
        delmem_complex_op()(temp);
    }
    delmem_complex_op()(hcc);
    delmem_complex_op()(scc);
    delmem_complex_op()(vcc);

    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diagH_subspace_init(hamilt::Hamilt<T, Device>* pHamilt,
    const T* psi,
    int psi_nr,
    int psi_nc,
    psi::Psi<T, Device>& evc,
    Real* en,
    const std::function<void(T*, const int)>& add_to_hcc,
    const std::function<void(const T* const, const int, const int)>& export_vcc)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_subspace_init");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace_init");

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
            "DiagoIterAssist::diagH_subspace_init",
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
            gemv_op<T, Device>()(ctx, 'C', psi_nc, nstart, &one, psi, psi_nc, hpsi, 1, &zero, hcc + i * nstart, 1);
        }

        T* spsi = temp;
        // do sPsi band by band
        for (int i = 0; i < nstart; i++)
        {
            syncmem_complex_op()(ppsi, psi + i * psi_nc, psi_nc);
            pHamilt->sPsi(ppsi, spsi, dmin, dmin, 1);

            gemv_op<T, Device>()(ctx,
                                 'C',
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

        gemm_op<T, Device>()(ctx, 'C', 'N', nstart, nstart, dmin, &one, ppsi, dmax, hpsi, dmax, &zero, hcc, nstart);

        T* spsi = temp;
        // do sPsi for all bands
        pHamilt->sPsi(ppsi, spsi, psi_temp.get_nbasis(), psi_temp.get_nbasis(), psi_temp.get_nbands());

        gemm_op<T, Device>()(ctx, 'C', 'N', nstart, nstart, dmin, &one, ppsi, dmax, spsi, dmax, &zero, scc, nstart);
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

    DiagoIterAssist::diagH_LAPACK(nstart, n_band, hcc, scc, nstart, en, vcc);

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
        gemm_op<T, Device>()(ctx,
                             'N',
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

        gemm_op<T, Device>()(ctx,
                             'N',
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

        // matrixSetToAnother<T, Device>()(ctx, n_band, evctemp, dmin, evc.get_pointer(), dmax);

        // delmem_complex_op()(ctx, evctemp);
    }

    delmem_complex_op()(hcc);
    delmem_complex_op()(scc);
    delmem_complex_op()(vcc);
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace_init");
}

template <typename T, typename Device>
void DiagoIterAssist<T, Device>::diagH_LAPACK(const int nstart,
                                              const int nbands,
                                              const T* hcc,
                                              const T* scc,
                                              const int ldh, // nstart
                                              Real* e,       // always in CPU
                                              T* vcc)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_LAPACK");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_LAPACK");

    Real* eigenvalues = nullptr;
    resmem_var_op()(eigenvalues, nstart);
    setmem_var_op()(eigenvalues, 0, nstart);

    dngvd_op<T, Device>()(ctx, nstart, ldh, hcc, scc, eigenvalues, vcc);

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

    ModuleBase::timer::tick("DiagoIterAssist", "diagH_LAPACK");
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

        T* hphi = temp;
        // do hPsi for all bands
        psi::Range all_bands_range(1, psi.get_current_k(), 0, nstart - 1);
        hpsi_info hpsi_in(&psi, all_bands_range, hphi);
        pHamilt->ops->hPsi(hpsi_in);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             &one,
                             psi.get_pointer(),
                             dmax,
                             hphi,
                             dmax,
                             &zero,
                             hcc,
                             nstart);

        T* sphi = temp;
        // do sPsi for all bands
        pHamilt->sPsi(psi.get_pointer(), sphi, dmax, dmin, nstart);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             &one,
                             psi.get_pointer(),
                             dmax,
                             sphi,
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
                                                const T* scc,
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
    DiagoIterAssist::diagH_LAPACK(nstart, nstart, hcc, scc, nstart, en, vcc);

    { // code block to calculate tar_mat
        gemm_op<T, Device>()(ctx,
                             'N',
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
                              const T* scc,
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
    DiagoIterAssist::diagH_LAPACK(nstart, nstart, hcc, scc, nstart, en, vcc);

    { // code block to calculate tar_mat
        const int dmin = evc.get_current_ngk();
        const int dmax = evc.get_nbasis();
        T* temp = nullptr;
        resmem_complex_op()(temp, nstart * dmax, "DiagSub::temp");
        setmem_complex_op()(temp, 0, nstart * dmax);
        gemm_op<T, Device>()(ctx,
                             'N',
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
        matrixSetToAnother<T, Device>()(ctx, n_band, temp, dmin, evc.get_pointer(), dmax);
        delmem_complex_op()(temp);
    }

    delmem_complex_op()(vcc);

    ModuleBase::timer::tick("DiagoIterAssist", "diag_subspace_psi");
}

template <typename T, typename Device>
bool DiagoIterAssist<T, Device>::test_exit_cond(const int& ntry, const int& notconv)
{
    //================================================================
    // If this logical function is true, need to do diagH_subspace
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
    // using diagH_subspace and cg method again. ntry++
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
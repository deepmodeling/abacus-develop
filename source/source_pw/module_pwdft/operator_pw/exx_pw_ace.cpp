#include "op_exx_pw.h"

namespace hamilt
{
template <typename T, typename Device>
void OperatorEXXPW<T, Device>::act_op_ace(const int nbands,
                                          const int nbasis,
                                          const int npol,
                                          const T *tmpsi_in,
                                          T *tmhpsi,
                                          const int ngk_ik,
                                          const bool is_first_node) const
{
    ModuleBase::timer::tick("OperatorEXXPW", "act_op_ace");
    //    std::cout << "act_op_ace" << std::endl;
    // hpsi += -Xi^\dagger * Xi * psi
    T* Xi_ace = Xi_ace_k[this->ik];
    int nbands_tot = psi.get_nbands();
    int nbasis_max = psi.get_nbasis();
    //    T* hpsi = nullptr;
    //    resmem_complex_op()(hpsi, nbands_tot * nbasis);
    //    setmem_complex_op()(hpsi, 0, nbands_tot * nbasis);
    T* Xi_psi = nullptr;
    resmem_complex_op()(Xi_psi, nbands_tot * nbands);
    setmem_complex_op()(Xi_psi, 0, nbands_tot * nbands);

    char trans_N = 'N', trans_T = 'T', trans_C = 'C';
    T intermediate_one = 1.0, intermediate_zero = 0.0, intermediate_minus_one = -1.0;
    // Xi * psi
    gemm_complex_op()(trans_N,
                      trans_N,
                      nbands_tot,
                      nbands,
                      nbasis,
                      &intermediate_one,
                      Xi_ace,
                      nbands_tot,
                      tmpsi_in,
                      nbasis,
                      &intermediate_zero,
                      Xi_psi,
                      nbands_tot
    );

    Parallel_Reduce::reduce_pool(Xi_psi, nbands_tot * nbands);

    // Xi^\dagger * (Xi * psi)
    gemm_complex_op()(trans_C,
                      trans_N,
                      nbasis,
                      nbands,
                      nbands_tot,
                      &intermediate_minus_one,
                      Xi_ace,
                      nbands_tot,
                      Xi_psi,
                      nbands_tot,
                      &intermediate_one,
                      tmhpsi,
                      nbasis
    );


    //    // negative sign, add to hpsi
    //    vec_add_vec_complex_op()(this->ctx, nbands * nbasis, tmhpsi, hpsi, -1, tmhpsi, 1);
    //    delmem_complex_op()(hpsi);
    delmem_complex_op()(Xi_psi);
    ModuleBase::timer::tick("OperatorEXXPW", "act_op_ace");

}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::construct_ace() const
{
    //    int nkb = p_exx_helper->psi.get_nbands() * p_exx_helper->psi.get_nk();
    int nbands = psi.get_nbands();
    int nbasis = psi.get_nbasis();
    int nk = psi.get_nk();

    int ik_save = this->ik;
    int * ik_ = const_cast<int*>(&this->ik);

    T intermediate_one = 1.0, intermediate_zero = 0.0;

    if (h_psi_ace == nullptr)
    {
        resmem_complex_op()(h_psi_ace, nbands * nbasis);
        setmem_complex_op()(h_psi_ace, 0, nbands * nbasis);
    }

    if (Xi_ace_k.size() != nk)
    {
        Xi_ace_k.resize(nk);
        for (int i = 0; i < nk; i++)
        {
            resmem_complex_op()(Xi_ace_k[i], nbands * nbasis);
        }
    }

    for (int i = 0; i < nk; i++)
    {
        setmem_complex_op()(Xi_ace_k[i], 0, nbands * nbasis);
    }

    if (L_ace == nullptr)
    {
        resmem_complex_op()(L_ace, nbands * nbands);
        setmem_complex_op()(L_ace, 0, nbands * nbands);
    }

    if (psi_h_psi_ace == nullptr)
    {
        resmem_complex_op()(psi_h_psi_ace, nbands * nbands);
    }

    if (first_iter) return;
    ModuleBase::timer::tick("OperatorEXXPW", "construct_ace");

    for (int ik = 0; ik < nk; ik++)
    {
        int npwk = wfcpw->npwk[ik];

        T* Xi_ace = Xi_ace_k[ik];
        psi.fix_kb(ik, 0);
        T* p_psi = psi.get_pointer();

        setmem_complex_op()(h_psi_ace, 0, nbands * nbasis);

        *ik_ = ik;

        act_op(
            nbands,
            nbasis,
            1,
            p_psi,
            h_psi_ace,
            nbasis,
            false
        );

        // psi_h_psi_ace = psi^\dagger * h_psi_ace
        // p_exx_helper->psi.fix_kb(0, 0);
        gemm_complex_op()('C',
                          'N',
                          nbands,
                          nbands,
                          npwk,
                          &intermediate_one,
                          p_psi,
                          nbasis,
                          h_psi_ace,
                          nbasis,
                          &intermediate_zero,
                          psi_h_psi_ace,
                          nbands);

        // reduction of psi_h_psi_ace, due to distributed memory
        Parallel_Reduce::reduce_pool(psi_h_psi_ace, nbands * nbands);

        T intermediate_minus_one = -1.0;
        axpy_complex_op()(nbands * nbands,
                          &intermediate_minus_one,
                          psi_h_psi_ace,
                          1,
                          L_ace,
                          1);


        int info = 0;
        char up = 'U', lo = 'L';

        lapack_potrf()(lo, nbands, L_ace, nbands);

        // expand for-loop
        for (int i = 0; i < nbands; ++i) {
            setmem_complex_op()(L_ace + i * nbands, 0, i);
        }

        // L_ace inv in place
        char non = 'N';
        lapack_trtri()(lo, non, nbands, L_ace, nbands);

        // Xi_ace = L_ace^-1 * h_psi_ace^dagger
        gemm_complex_op()('N',
                          'C',
                          nbands,
                          npwk,
                          nbands,
                          &intermediate_one,
                          L_ace,
                          nbands,
                          h_psi_ace,
                          nbasis,
                          &intermediate_zero,
                          Xi_ace,
                          nbands);

        // clear mem
        setmem_complex_op()(h_psi_ace, 0, nbands * nbasis);
        setmem_complex_op()(psi_h_psi_ace, 0, nbands * nbands);
        setmem_complex_op()(L_ace, 0, nbands * nbands);

    }

    *ik_ = ik_save;
    ModuleBase::timer::tick("OperatorEXXPW", "construct_ace");

}

template <typename T, typename Device>
double OperatorEXXPW<T, Device>::cal_exx_energy_ace(psi::Psi<T, Device>* ppsi_) const
{
    double Eexx = 0;

    psi::Psi<T, Device> psi_ = *ppsi_;
    int* ik_ = const_cast<int*>(&this->ik);
    int ik_save = this->ik;
    for (int i = 0; i < wfcpw->nks; i++)
    {
        setmem_complex_op()(h_psi_ace, 0, psi_.get_nbands() * psi_.get_nbasis());
        *ik_ = i;
        psi_.fix_kb(i, 0);
        T* psi_i = psi_.get_pointer();
        act_op_ace(psi_.get_nbands(), psi_.get_nbasis(), 1, psi_i, h_psi_ace, 0, true);

        for (int nband = 0; nband < psi_.get_nbands(); nband++)
        {
            psi_.fix_kb(i, nband);
            T* psi_i_n = psi_.get_pointer();
            T* hpsi_i_n = h_psi_ace + nband * psi_.get_nbasis();
            double wg_i_n = (*wg)(i, nband);
            // Eexx += dot(psi_i_n, h_psi_i_n)
            Eexx += dot_op()(psi_.get_nbasis(), psi_i_n, hpsi_i_n, false) * wg_i_n * 2;
        }
    }

    Parallel_Reduce::reduce_pool(Eexx);
    *ik_ = ik_save;
    return Eexx;
}
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
}
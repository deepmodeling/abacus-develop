#include "op_pw_exx.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_xc/exx_info.h"

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
    ModuleBase::timer::start("OperatorEXXPW", "act_op_ace");
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

#ifdef __MPI
    Parallel_Common::reduce_dev<T, Device>(Xi_psi, nbands_tot * nbands, POOL_WORLD);
#endif

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

    delmem_complex_op()(Xi_psi);
    ModuleBase::timer::end("OperatorEXXPW", "act_op_ace");

}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::construct_ace() const
{
    int nbands = psi.get_nbands();
    int nbasis = psi.get_nbasis();
    int nk = psi.get_nk();
    const int hpsi_size = nbands * nbasis;

    int* ik_ = const_cast<int*>(&this->ik);
    int ik_save = this->ik;

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
    ModuleBase::timer::start("OperatorEXXPW", "construct_ace");

    int nk_max = kv->para_k.get_max_nks_pool();
    int nspin_fac = nspin_ == 2 ? 2 : 1;
    for (int ispin = 0; ispin < nspin_fac; ispin++)
    {
        for (int ik0 = 0; ik0 < nk_max; ik0++)
        {
            const int nks_per_spin = wfcpw->nks / nspin_fac;
            const bool has_kpoint = ik0 < nks_per_spin;
            // Use an out-of-range sentinel on idle pools. They still enter the
            // same q-point collectives, but never access local k-point data.
            const int ik = has_kpoint ? ik0 + ispin * nks_per_spin : wfcpw->nks;
            *ik_ = ik;
            T* p_psi = nullptr;
            if (has_kpoint)
            {
                psi.fix_kb(ik, 0);
                p_psi = psi.get_pointer();
            }
            setmem_complex_op()(h_psi_ace, 0, hpsi_size);
            act_op_kpar(nbands, nbasis, p_psi, h_psi_ace, ispin, has_kpoint);

            if (has_kpoint)
            {
                const int npwk = wfcpw->npwk[ik];
                T* Xi_ace = Xi_ace_k[ik];
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
#ifdef __MPI
                Parallel_Common::reduce_dev<T, Device>(psi_h_psi_ace, nbands * nbands, POOL_WORLD);
#endif

                T intermediate_minus_one = -1.0;
                axpy_complex_op()(nbands * nbands,
                                  &intermediate_minus_one,
                                  psi_h_psi_ace,
                                  1,
                                  L_ace,
                                  1);


                int info = 0;
                char up = 'U', lo = 'L';

                // for (int i = 0; i < nbands; ++i)
                // {
                //     for (int j = 0; j < nbands; ++j)
                //     {
                //         // std::cout << L_ace[i * nbands + j]. << " ";
                //         if (L_ace[i * nbands + j].imag() >= 0.0)
                //         {
                //             std::cout << L_ace[i * nbands + j].real() << "+" << L_ace[i * nbands + j].imag() << "im ";
                //         }
                //         else
                //         {
                //             std::cout << L_ace[i * nbands + j].real() << L_ace[i * nbands + j].imag() << "im ";
                //         }
                //     }
                //     std::cout << ";" << std::endl;
                // }
                // MPI_Barrier(MPI_COMM_WORLD);
                // MPI_Abort(MPI_COMM_WORLD, 0);

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
                setmem_complex_op()(h_psi_ace, 0, hpsi_size);
                setmem_complex_op()(psi_h_psi_ace, 0, nbands * nbands);
                setmem_complex_op()(L_ace, 0, nbands * nbands);
            }
        }
    }

    *ik_ = ik_save;

    ModuleBase::timer::end("OperatorEXXPW", "construct_ace");

}

template <typename T, typename Device>
double OperatorEXXPW<T, Device>::cal_exx_energy_ace(psi::Psi<T, Device>* ppsi_) const
{
    double Eexx = 0;
    int nspin_fac = nspin_ == 2 ? 2 : 1;
    psi::Psi<T, Device> psi_ = *ppsi_;
    int* ik_ = const_cast<int*>(&this->ik);
    int ik_save = this->ik;
    Real hybrid_alpha = this->hybrid_alpha;
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
            Eexx += dot_op()(psi_.get_nbasis(), psi_i_n, hpsi_i_n, false) * wg_i_n;
        }
    }

    Parallel_Reduce::reduce_all(Eexx);
    *ik_ = ik_save;
    Eexx = Eexx / hybrid_alpha / 2; // This factor of 2 is from the definition of EXX energy.
    return Eexx;
}

// Explicit instantiation for members defined in this translation unit.
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt

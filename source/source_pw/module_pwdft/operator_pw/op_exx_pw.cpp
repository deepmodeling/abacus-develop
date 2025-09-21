#include "op_exx_pw.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/klist.h"
#include "source_hamilt/operator.h"
#include "source_psi/psi.h"
#include "source_pw/module_pwdft/global.h"
#include "source_pw/module_pwdft/kernels/cal_density_real_op.h"
#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_pw/module_pwdft/kernels/mul_potential_op.h"
#include "source_pw/module_pwdft/kernels/vec_mul_vec_complex_op.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>

namespace hamilt
{

template <typename T, typename Device>
OperatorEXXPW<T, Device>::OperatorEXXPW(const int* isk_in,
                                        const ModulePW::PW_Basis_K* wfcpw_in,
                                        const ModulePW::PW_Basis* rhopw_in,
                                        K_Vectors *kv_in,
                                        const UnitCell *ucell)
    : isk(isk_in), wfcpw(wfcpw_in), rhopw(rhopw_in), kv(kv_in), ucell(ucell)
{
    gamma_extrapolation = PARAM.inp.exx_gamma_extrapolation;
    bool is_mp = kv_in->get_is_mp();
#ifdef __MPI
    Parallel_Common::bcast_bool(is_mp);
#endif
    if (!is_mp)
    {
        gamma_extrapolation = false;
    }
    if (GlobalV::KPAR != 1)
    {
        // GlobalV::ofs_running << "EXX Calculation does not support k-point parallelism" << std::endl;
        ModuleBase::WARNING_QUIT("OperatorEXXPW", "EXX Calculation does not support k-point parallelism");
    }

    this->classname = "OperatorEXXPW";
    this->ctx = nullptr;
    this->cpu_ctx = nullptr;
    this->cal_type = hamilt::calculation_type::pw_exx;

    // allocate real space memory
    // assert(wfcpw->nrxx == rhopw->nrxx);
    resmem_complex_op()(psi_nk_real, wfcpw->nrxx);
    resmem_complex_op()(psi_mq_real, wfcpw->nrxx);
    resmem_complex_op()(density_real, rhopw->nrxx);
    resmem_complex_op()(h_psi_real, rhopw->nrxx);
    // allocate density recip space memory
    resmem_complex_op()(density_recip, rhopw->npw);
    // allocate h_psi recip space memory
    resmem_complex_op()(h_psi_recip, wfcpw->npwk_max);
    // resmem_complex_op()(this->ctx, psi_all_real, wfcpw->nrxx * GlobalV::NBANDS);

    int nks = wfcpw->nks;
    int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    resmem_real_op()(pot, rhopw->npw * nks * nks);

    tpiba = ucell->tpiba;
    Real tpiba2 = tpiba * tpiba;

    rhopw_dev = new ModulePW::PW_Basis(wfcpw->get_device(), rhopw->get_precision());
    rhopw_dev->fft_bundle.setfft(wfcpw->get_device(), rhopw->get_precision());
#ifdef __MPI
    rhopw_dev->initmpi(rhopw->poolnproc, rhopw->poolrank, rhopw->pool_world);
#endif
    // here we can actually use different ecut to init the grids
    rhopw_dev->initgrids(rhopw->lat0, rhopw->latvec, rhopw->gridecut_lat * rhopw->tpiba2);
    rhopw_dev->initgrids(rhopw->lat0, rhopw->latvec, rhopw->nx, rhopw->ny, rhopw->nz);
    rhopw_dev->initparameters(rhopw->gamma_only, rhopw->ggecut * rhopw->tpiba2, rhopw->distribution_type, rhopw->xprime);
    rhopw_dev->setuptransform();
    rhopw_dev->collect_local_pw();

}   // end of constructor

template <typename T, typename Device>
OperatorEXXPW<T, Device>::~OperatorEXXPW()
{
    // use delete_memory_op to delete the allocated pws
    delmem_complex_op()(psi_nk_real);
    delmem_complex_op()(psi_mq_real);
    delmem_complex_op()(density_real);
    delmem_complex_op()(h_psi_real);
    delmem_complex_op()(density_recip);
    delmem_complex_op()(h_psi_recip);

    delmem_real_op()(pot);

    delmem_complex_op()(h_psi_ace);
    delmem_complex_op()(psi_h_psi_ace);
    delmem_complex_op()(L_ace);
    for (auto &Xi_ace: Xi_ace_k)
    {
        delmem_complex_op()(Xi_ace);
    }
    Xi_ace_k.clear();
    delete rhopw_dev;
}

template <typename T>
inline bool is_finite(const T &val)
{
    return std::isfinite(val);
}

template <>
inline bool is_finite(const std::complex<float> &val)
{
    return std::isfinite(val.real()) && std::isfinite(val.imag());
}

template <>
inline bool is_finite(const std::complex<double> &val)
{
    return std::isfinite(val.real()) && std::isfinite(val.imag());
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::act(const int nbands,
                                   const int nbasis,
                                   const int npol,
                                   const T *tmpsi_in,
                                   T *tmhpsi,
                                   const int ngk_ik,
                                   const bool is_first_node) const
{
    if (first_iter) return;

    if (is_first_node)
    {
        setmem_complex_op()(tmhpsi, 0, nbasis*nbands/npol);
    }

    if (PARAM.inp.exxace && GlobalC::exx_info.info_global.separate_loop)
    {
        act_op_ace(nbands, nbasis, npol, tmpsi_in, tmhpsi, ngk_ik, is_first_node);
    }
    else
    {
        act_op(nbands, nbasis, npol, tmpsi_in, tmhpsi, ngk_ik, is_first_node);
    }
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::act_op(const int nbands,
                                   const int nbasis,
                                   const int npol,
                                   const T *tmpsi_in,
                                   T *tmhpsi,
                                   const int ngk_ik,
                                   const bool is_first_node) const
{
//    std::cout << "nbands: " << nbands
//              << " nbasis: " << nbasis
//              << " npol: " << npol
//              << " ngk_ik: " << ngk_ik
//              << " is_first_node: " << is_first_node
//              << std::endl;
    if (!potential_got)
    {
        get_exx_potential<Real, Device>(kv, wfcpw, rhopw_dev, pot, tpiba, gamma_extrapolation, ucell->omega);
        potential_got = true;
    }

//    set_psi(&p_exx_helper->psi);

    ModuleBase::timer::tick("OperatorEXXPW", "act_op");

    setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
    setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_recip, 0, rhopw_dev->npw);
    // setmem_complex_op()(psi_all_real, 0, wfcpw->nrxx * GlobalV::NBANDS);
    // std::map<std::pair<int, int>, bool> has_real;
    setmem_complex_op()(psi_nk_real, 0, wfcpw->nrxx);
    setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);

    // ik fixed here, select band n
    for (int n_iband = 0; n_iband < nbands; n_iband++)
    {
        const T *psi_nk = tmpsi_in + n_iband * nbasis;
        // retrieve \psi_nk in real space
        wfcpw->recip_to_real(ctx, psi_nk, psi_nk_real, this->ik);

        // for \psi_nk, get the pw of iq and band m
        auto q_points = get_q_points(this->ik);
        Real nqs = q_points.size();
        for (int iq: q_points)
        {
//            std::cout << "ik" << this->ik << " iq" << iq << std::endl;
            for (int m_iband = 0; m_iband < psi.get_nbands(); m_iband++)
            {
                // double wg_mqb_real = GlobalC::exx_helper.wg(iq, m_iband);
                double wg_mqb_real = (*wg)(this->ik, m_iband);
                T wg_mqb = wg_mqb_real;
                if (wg_mqb_real < 1e-12)
                {
                    continue;
                }

                const T* psi_mq = get_pw(m_iband, iq);
                wfcpw->recip_to_real(ctx, psi_mq, psi_mq_real, iq);

                // direct multiplication in real space, \psi_nk(r) * \psi_mq(r)
                cal_density_recip(psi_nk_real, psi_mq_real, ucell->omega);

                // bring the density to recip space
                // rhopw->real2recip(density_real, density_recip);

                // multiply the density with the potential in recip space
                multiply_potential(density_recip, this->ik, iq);

                // bring the potential back to real space
                // rhopw_dev->recip2real(density_recip, density_real);
                rho_recip2real(density_recip, density_real);

                if (false)
                {
                    // do nothing
                }
                else
                {
                    vec_mul_vec_complex_op<T, Device>()(density_real, psi_mq_real, density_real, wfcpw->nrxx);
                }

                T wk_iq = kv->wk[iq];
                T wk_ik = kv->wk[this->ik];

                T tmp_scalar = wg_mqb / wk_iq / nqs;
                axpy_complex_op()(wfcpw->nrxx,
                                  &tmp_scalar,
                                  density_real,
                                  1,
                                  h_psi_real,
                                  1);

            } // end of m_iband
            setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
            setmem_complex_op()(density_recip, 0, rhopw_dev->npw);
            setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);

        } // end of iq
        T* h_psi_nk = tmhpsi + n_iband * nbasis;
        Real hybrid_alpha = GlobalC::exx_info.info_global.hybrid_alpha;
        wfcpw->real_to_recip(ctx, h_psi_real, h_psi_nk, this->ik, true, hybrid_alpha);
        setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);

    }

    ModuleBase::timer::tick("OperatorEXXPW", "act_op");

}

template <typename T, typename Device>
std::vector<int> OperatorEXXPW<T, Device>::get_q_points(const int ik) const
{
    // stored in q_points
    if (q_points.find(ik) != q_points.end())
    {
        return q_points.find(ik)->second;
    }

    std::vector<int> q_points_ik;

    // if () // downsampling
    {
        for (int iq = 0; iq < wfcpw->nks; iq++)
        {
            if (PARAM.inp.nspin ==1 )
            {
                q_points_ik.push_back(iq);
            }
            else if (PARAM.inp.nspin == 2)
            {
                int nk_fac = 2;
                int nk = wfcpw->nks / nk_fac;
                if (iq / nk == ik / nk)
                {
                    q_points_ik.push_back(iq);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("OperatorEXXPW", "nspin == 4 not supported");
            }
        }
    }
    // else
    // {
    //     for (int iq = 0; iq < wfcpw->nks; iq++)
    //     {
    //         kv->
    //     }
    // }

    q_points[ik] = q_points_ik;
    return q_points_ik;
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::multiply_potential(T *density_recip, int ik, int iq) const
{
    ModuleBase::timer::tick("OperatorEXXPW", "multiply_potential");
    int npw = rhopw_dev->npw;
    int nks = wfcpw->nks;
    int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    int nk = nks / nk_fac;

    mul_potential_op<T, Device>()(pot, density_recip, npw, nks, ik, iq);

    ModuleBase::timer::tick("OperatorEXXPW", "multiply_potential");
}

template <typename T, typename Device>
const T *OperatorEXXPW<T, Device>::get_pw(const int m, const int iq) const
{
    // return pws[iq].get() + m * wfcpw->npwk[iq];
    psi.fix_kb(iq, m);
    T* psi_mq = psi.get_pointer();
    return psi_mq;
}

template <typename T, typename Device>
template <typename T_in, typename Device_in>
OperatorEXXPW<T, Device>::OperatorEXXPW(const OperatorEXXPW<T_in, Device_in> *op)
{
    // copy all the datas
    this->isk = op->isk;
    this->wfcpw = op->wfcpw;
    this->rhopw = op->rhopw;
    this->rhopw_dev = op->rhopw_dev;
    this->psi = op->psi;
    this->ctx = op->ctx;
    this->cpu_ctx = op->cpu_ctx;
    resmem_complex_op()(this->ctx, psi_nk_real, wfcpw->nrxx);
    resmem_complex_op()(this->ctx, psi_mq_real, wfcpw->nrxx);
    resmem_complex_op()(this->ctx, density_real, rhopw_dev->nrxx);
    resmem_complex_op()(this->ctx, h_psi_real, rhopw_dev->nrxx);
    resmem_complex_op()(this->ctx, density_recip, rhopw_dev->npw);
    resmem_complex_op()(this->ctx, h_psi_recip, wfcpw->npwk_max);
//    this->pws.resize(wfcpw->nks);


}

template <typename T, typename Device>
double OperatorEXXPW<T, Device>::cal_exx_energy(psi::Psi<T, Device> *psi_) const
{
    if (PARAM.inp.exxace && GlobalC::exx_info.info_global.separate_loop)
    {
        return cal_exx_energy_ace(psi_);
    }
    else
    {
        return cal_exx_energy_op(psi_);
    }
}

template <typename T, typename Device>
double OperatorEXXPW<T, Device>::cal_exx_energy_op(psi::Psi<T, Device> *ppsi_) const
{
    psi::Psi<T, Device> psi_ = *ppsi_;

    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    setmem_complex_op()(psi_nk_real, 0, wfcpw->nrxx);
    setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);
    setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
    setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_recip, 0, rhopw_dev->npw);

    if (wg == nullptr) return 0.0;
    const int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    double Eexx_ik_real = 0.0;
    for (int ik = 0; ik < wfcpw->nks; ik++)
    {
        //        auto k = this->pw_wfc->kvec_c[ik];
        //        std::cout << k << std::endl;
        for (int n_iband = 0; n_iband < psi.get_nbands(); n_iband++)
        {
            setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
            setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
            setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
            setmem_complex_op()(density_recip, 0, rhopw_dev->npw);

            // double wg_ikb_real = GlobalC::exx_helper.wg(this->ik, n_iband);
            double wg_ikb_real = (*wg)(ik, n_iband);
            T wg_ikb = wg_ikb_real;
            if (wg_ikb_real < 1e-12)
            {
                continue;
            }

            // const T *psi_nk = get_pw(n_iband, ik);
            psi.fix_kb(ik, n_iband);
            const T* psi_nk = psi.get_pointer();
            // retrieve \psi_nk in real space
            wfcpw->recip_to_real(ctx, psi_nk, psi_nk_real, ik);

            // for \psi_nk, get the pw of iq and band m
            // q_points is a vector of integers, 0 to nks-1
            std::vector<int> q_points;
            for (int iq = 0; iq < wfcpw->nks; iq++)
            {
                q_points.push_back(iq);
            }
            double nqs = q_points.size();

            for (int iq: q_points)
            {
                for (int m_iband = 0; m_iband < psi.get_nbands(); m_iband++)
                {
                    // double wg_f = GlobalC::exx_helper.wg(iq, m_iband);
                    double wg_iqb_real = (*wg)(iq, m_iband);
                    T wg_iqb = wg_iqb_real;
                    if (wg_iqb_real < 1e-12)
                    {
                        continue;
                    }

                    psi_.fix_kb(iq, m_iband);
                    const T* psi_mq = psi_.get_pointer();
                    // const T* psi_mq = get_pw(m_iband, iq);
                    wfcpw->recip_to_real(ctx, psi_mq, psi_mq_real, iq);

                    cal_density_recip(psi_nk_real, psi_mq_real, ucell->omega);

                    int nks = wfcpw->nks;
                    int npw = rhopw_dev->npw;
                    int nk = nks / nk_fac;
                    Eexx_ik_real += exx_cal_energy_op<T, Device>()(density_recip, pot + ik * nks * npw + iq * npw, wg_iqb_real / nqs * wg_ikb_real / kv->wk[ik], npw);

                } // m_iband

            } // iq

        } // n_iband

    } // ik
    Eexx_ik_real *= 0.5 * ucell->omega;
    Parallel_Reduce::reduce_pool(Eexx_ik_real);
    //    std::cout << "omega = " << this_->pelec->omega << " tpiba = " << this_->pw_rho->tpiba2 << " exx_div = " << exx_div << std::endl;

    setmem_complex_op()(psi_nk_real, 0, wfcpw->nrxx);
    setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);
    setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
    setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_recip, 0, rhopw_dev->npw);

    return Eexx_ik_real;
}

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>::cal_density_recip(const std::complex<double>* psi_nk_real,
                                                                                const std::complex<double>* psi_mq_real,
                                                                                double omega) const
{
    cal_density_real_op<std::complex<double>, base_device::DEVICE_CPU>()(psi_nk_real, psi_mq_real, density_real, omega, wfcpw->nrxx);
    rhopw_dev->real2recip(density_real, density_recip);
}

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>::cal_density_recip(const std::complex<float>* psi_nk_real,
                                                                                const std::complex<float>* psi_mq_real,
                                                                                double omega) const
{
    cal_density_real_op<std::complex<float>, base_device::DEVICE_CPU>()(psi_nk_real, psi_mq_real, density_real, omega, wfcpw->nrxx);
    rhopw_dev->real2recip(density_real, density_recip);
}

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>::rho_recip2real(const std::complex<double>* rho_recip,
                                                                             std::complex<double>* rho_real,
                                                                             bool add,
                                                                             double factor) const
{
    rhopw_dev->recip2real(rho_recip, rho_real, add, factor);
}

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>::rho_recip2real(const std::complex<float>* rho_recip,
                                                                             std::complex<float>* rho_real,
                                                                             bool add,
                                                                             float factor) const
{
    rhopw_dev->recip2real(rho_recip, rho_real, add, factor);
}

template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>;
template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>;

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>::cal_density_recip(const std::complex<double>* psi_nk_real,
                                                                                const std::complex<double>* psi_mq_real,
                                                                                double omega) const
{
    cal_density_real_op<std::complex<double>, base_device::DEVICE_GPU>()(psi_nk_real, psi_mq_real, density_real, omega, wfcpw->nrxx);
    rhopw_dev->real2recip_gpu(density_real, density_recip);
}

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>::cal_density_recip(const std::complex<float>* psi_nk_real,
                                                                                const std::complex<float>* psi_mq_real,
                                                                                double omega) const
{
    cal_density_real_op<std::complex<float>, base_device::DEVICE_GPU>()(psi_nk_real, psi_mq_real, density_real, omega, wfcpw->nrxx);
    rhopw_dev->real2recip_gpu(density_real, density_recip);
}

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>::rho_recip2real(const std::complex<double>* rho_recip,
                                                                             std::complex<double>* rho_real,
                                                                             bool add,
                                                                             double factor) const
{
    rhopw_dev->recip2real_gpu(rho_recip, rho_real, add, factor);
}

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>::rho_recip2real(const std::complex<float>* rho_recip,
                                                                             std::complex<float>* rho_real,
                                                                             bool add,
                                                                             float factor) const
{
    rhopw_dev->recip2real_gpu(rho_recip, rho_real, add, factor);
}

#endif

} // namespace hamilt

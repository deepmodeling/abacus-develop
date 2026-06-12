#include "source_hamilt/module_xc/exx_info.h"
#include "op_pw_exx.h"
#include "source_pw/module_pwdft/kernels/cal_density_real_op.h"
#include "source_pw/module_pwdft/kernels/exx_q_state_op.h"
#include "source_pw/module_pwdft/kernels/exx_stress_op.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_reduce.h"
#include "stress_pw.h"

#include <cmath>
#include <algorithm>
#include <type_traits>
#include <vector>

#if defined(__ROCM) && !defined(__CUDA)
#error "PW EXX GPU stress implementation is not implemented for ROCm in this merge."
#endif

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_exx(ModuleBase::matrix& sigma,
                                           const ModuleBase::matrix& wg,
                                           ModulePW::PW_Basis* rhopw,
                                           ModulePW::PW_Basis_K* wfcpw,
                                           const K_Vectors *p_kv,
                                           const psi::Psi <std::complex<FPTYPE>, Device>* d_psi_in, const UnitCell& ucell)
{
    bool gamma_extrapolation = PARAM.inp.exx_gamma_extrapolation;
    bool is_mp = p_kv->get_is_mp();
#ifdef __MPI
    Parallel_Common::bcast_bool(is_mp);
#endif
    if (!is_mp)
    {
        gamma_extrapolation = false;
    }

    // T is complex of FPTYPE, if FPTYPE is double, T is std::complex<double>
    // but if FPTYPE is std::complex<double>, T is still std::complex<double>
    using T = std::complex<FPTYPE>;
    using Real = FPTYPE;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using resmem_real_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using delmem_real_op = base_device::memory::delete_memory_op<Real, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_real_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_real_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;

    if (GlobalV::KPAR != 1 && !(PARAM.inp.exxace && GlobalC::exx_info.info_global.separate_loop))
    {
        ModuleBase::WARNING_QUIT("Stress_PW::stress_exx",
                                 "PW EXX KPAR stress is supported only with exxace=1 and exx_separate_loop=1");
    }

    const int nspin_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    const Real k_spin_degeneracy = PARAM.inp.nspin == 1 ? 2.0 : 1.0;
    double omega = ucell.omega;
    double tpiba = ucell.tpiba;
    double tpiba2 = ucell.tpiba2;
    double omega_inv = 1.0 / omega;

    // allocate space
    T* psi_nk_real = nullptr;
    T* psi_mq_real = nullptr;
    T* density_real = nullptr;
    T* density_recip = nullptr;
    Real* pot = nullptr; // This factor is 2x of the potential in 10.1103/PhysRevB.73.125120
    Real* pot_stress = nullptr;
    Real* sigma_exx_device = nullptr;
    Real* gcar_flat = nullptr;
    ModulePW::PW_Basis* rhopw_exx_owned = nullptr;
    ModulePW::PW_Basis* rhopw_exx = rhopw;
    const bool use_device_q_tile = PARAM.inp.exx_use_q_tile
                                   && !std::is_same<Device, base_device::DEVICE_CPU>::value;
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value && !use_device_q_tile)
    {
        ModuleBase::WARNING_QUIT("Stress_PW::stress_exx",
                                 "GPU PW EXX stress requires exx_use_q_tile=1");
    }
#if defined(__ROCM) && !defined(__CUDA)
    if (use_device_q_tile)
    {
        ModuleBase::WARNING_QUIT("Stress_PW::stress_exx",
                                 "GPU q-tile PW EXX stress path is not implemented for ROCm yet");
    }
#endif
    const bool use_gpu_q_tile = use_device_q_tile;
    if (use_gpu_q_tile)
    {
        double ecut_exx = PARAM.inp.ecutexx;
        if (ecut_exx == 0.0)
        {
            ecut_exx = PARAM.inp.ecutrho;
        }
        rhopw_exx_owned = new ModulePW::PW_Basis(wfcpw->get_device(), rhopw->get_precision());
        rhopw_exx_owned->fft_bundle.setfft(wfcpw->get_device(), rhopw->get_precision());
#ifdef __MPI
        rhopw_exx_owned->initmpi(rhopw->poolnproc, rhopw->poolrank, rhopw->pool_world);
#endif
        rhopw_exx_owned->initgrids(rhopw->lat0, rhopw->latvec, ecut_exx);
        rhopw_exx_owned->initgrids(rhopw->lat0, rhopw->latvec, rhopw->nx, rhopw->ny, rhopw->nz);
        rhopw_exx_owned->initparameters(rhopw->gamma_only, ecut_exx, rhopw->distribution_type, rhopw->xprime);
        rhopw_exx_owned->setuptransform();
        rhopw_exx_owned->collect_local_pw();
        rhopw_exx = rhopw_exx_owned;
    }

    resmem_complex_op()(psi_nk_real, wfcpw->nrxx);
    resmem_complex_op()(psi_mq_real, wfcpw->nrxx);
    resmem_complex_op()(density_real, rhopw_exx->nrxx);
    resmem_complex_op()(density_recip, rhopw_exx->npw);
    resmem_real_op()(pot, rhopw_exx->npw);
    resmem_real_op()(pot_stress, rhopw_exx->npw);
    resmem_real_op()(sigma_exx_device, 6);

    if (use_gpu_q_tile)
    {
        std::vector<Real> gcar_host(static_cast<std::size_t>(rhopw_exx->npw) * 3);
        for (int ig = 0; ig < rhopw_exx->npw; ++ig)
        {
            gcar_host[ig * 3] = static_cast<Real>(rhopw_exx->gcar[ig].x);
            gcar_host[ig * 3 + 1] = static_cast<Real>(rhopw_exx->gcar[ig].y);
            gcar_host[ig * 3 + 2] = static_cast<Real>(rhopw_exx->gcar[ig].z);
        }
        resmem_real_op()(gcar_flat, gcar_host.size());
        syncmem_real_h2d_op()(gcar_flat, gcar_host.data(), gcar_host.size());
        base_device::memory::set_memory_op<Real, Device>()(sigma_exx_device, 0, 6);
    }

    // hamilt::get_exx_potential<Real, Device>(p_kv, wfcpw, rhopw, pot, tpiba, gamma_extrapolation, omega);
    // hamilt::get_exx_stress_potential<Real, Device>(p_kv, wfcpw, rhopw, pot_stress, tpiba, gamma_extrapolation, omega);

    auto load_q_real = [&](const K_Vectors::ExxFullQPoint& qpoint, int ispin_in, int mband, T* out) {
        const int iq_rep_spin = p_kv->exx_rep_spin_index(qpoint, ispin_in);
        d_psi_in->fix_kb(iq_rep_spin, mband);
        T* psi_mq = d_psi_in->get_pointer();
        if (qpoint.identity || qpoint.conjugate_only)
        {
            wfcpw->template recip_to_real<T, Device>(psi_mq, out, iq_rep_spin);
            if (qpoint.conjugate_only)
            {
                hamilt::exx_conjugate_real_op<T, Device>()(out, out, wfcpw->nrxx);
            }
        }
        else
        {
            if (std::is_same<Device, base_device::DEVICE_CPU>::value)
            {
                const auto remap = hamilt::build_exx_symmetry_remap(wfcpw, qpoint, iq_rep_spin, false);
                if (qpoint.time_reversal)
                {
                    wfcpw->recip2real_remapped_conjugate(psi_mq,
                                                         out,
                                                         static_cast<int>(remap.rep_igl.size()),
                                                         remap.rep_igl.data(),
                                                         remap.fft_isz.data(),
                                                         remap.phase.data());
                }
                else
                {
                    wfcpw->recip2real_remapped(psi_mq,
                                               out,
                                               static_cast<int>(remap.rep_igl.size()),
                                               remap.rep_igl.data(),
                                               remap.fft_isz.data(),
                                               remap.phase.data());
                }
            }
            else
            {
                const auto remap = hamilt::build_exx_symmetry_remap(wfcpw, qpoint, iq_rep_spin, true);
                if (qpoint.time_reversal)
                {
                    wfcpw->recip2real_remapped_conjugate(psi_mq,
                                                         out,
                                                         static_cast<int>(remap.rep_igl.size()),
                                                         remap.rep_igl.data(),
                                                         remap.fft_isz.data(),
                                                         remap.phase.data());
                }
                else
                {
                    wfcpw->recip2real_remapped(psi_mq,
                                               out,
                                               static_cast<int>(remap.rep_igl.size()),
                                               remap.rep_igl.data(),
                                               remap.fft_isz.data(),
                                               remap.phase.data());
                }
            }
        }
    };

    const std::vector<const K_Vectors::ExxFullQPoint*> q_points = [&]() {
        std::vector<const K_Vectors::ExxFullQPoint*> points;
        for (const auto& qpoint: p_kv->exx_full_q_map)
        {
            if (qpoint.active)
            {
                points.push_back(&qpoint);
            }
        }
        return points;
    }();

    const int nbands_psi = d_psi_in->get_nbands();
    const int source_tile_size = std::max(1, std::min(PARAM.inp.exx_band_tile_size, nbands_psi));
    const int q_tile_size = q_points.empty() ? 1
                                             : std::max(1, std::min(PARAM.inp.exx_q_tile_size,
                                                                    static_cast<int>(q_points.size())));
    T* q_real_tile = nullptr;
    std::vector<Real> q_weights;
    if (PARAM.inp.exx_use_q_tile)
    {
        resmem_complex_op()(q_real_tile,
                            static_cast<std::size_t>(q_tile_size) * static_cast<std::size_t>(source_tile_size)
                                * static_cast<std::size_t>(wfcpw->nrxx));
        q_weights.resize(static_cast<std::size_t>(q_tile_size) * static_cast<std::size_t>(source_tile_size), 0);
    }

    for (int ispin = 0; ispin < nspin_fac; ++ispin)
    {
        for (const auto& kpoint: p_kv->exx_full_k_map)
        {
            if (!kpoint.active)
            {
                continue;
            }
            const int ik_rep_spin = p_kv->exx_rep_spin_index(kpoint, ispin);
            const bool own_kpoint = kpoint.rep_pool == GlobalV::MY_POOL;
            for (int nband = 0; nband < d_psi_in->get_nbands(); nband++)
            {
                double wg_nkb = 0.0;
                double wk_nk = 0.0;
                if (own_kpoint)
                {
                    wg_nkb = wg(ik_rep_spin, nband);
                    wk_nk = p_kv->wk[ik_rep_spin];
                }
#ifdef __MPI
                MPI_Bcast(&wg_nkb, 1, MPI_DOUBLE, p_kv->para_k.get_startpro_pool(kpoint.rep_pool), MPI_COMM_WORLD);
                MPI_Bcast(&wk_nk, 1, MPI_DOUBLE, p_kv->para_k.get_startpro_pool(kpoint.rep_pool), MPI_COMM_WORLD);
#endif
                if (wg_nkb >= 1e-12 && own_kpoint)
                {
                    d_psi_in->fix_kb(ik_rep_spin, nband);
                    T* psi_nk = d_psi_in->get_pointer();
                    if (kpoint.identity || kpoint.conjugate_only)
                    {
                        wfcpw->template recip_to_real<T, Device>(psi_nk, psi_nk_real, ik_rep_spin);
                        if (kpoint.conjugate_only)
                        {
                            hamilt::exx_conjugate_real_op<T, Device>()(psi_nk_real, psi_nk_real, wfcpw->nrxx);
                        }
                    }
                    else
                    {
                        if (std::is_same<Device, base_device::DEVICE_CPU>::value)
                        {
                            const auto remap = hamilt::build_exx_symmetry_remap(wfcpw, kpoint, ik_rep_spin, false);
                            if (kpoint.time_reversal)
                            {
                                wfcpw->recip2real_remapped_conjugate(psi_nk,
                                                                     psi_nk_real,
                                                                     static_cast<int>(remap.rep_igl.size()),
                                                                     remap.rep_igl.data(),
                                                                     remap.fft_isz.data(),
                                                                     remap.phase.data());
                            }
                            else
                            {
                                wfcpw->recip2real_remapped(psi_nk,
                                                           psi_nk_real,
                                                           static_cast<int>(remap.rep_igl.size()),
                                                           remap.rep_igl.data(),
                                                           remap.fft_isz.data(),
                                                           remap.phase.data());
                            }
                        }
                        else
                        {
                            const auto remap = hamilt::build_exx_symmetry_remap(wfcpw, kpoint, ik_rep_spin, true);
                            if (kpoint.time_reversal)
                            {
                                wfcpw->recip2real_remapped_conjugate(psi_nk,
                                                                     psi_nk_real,
                                                                     static_cast<int>(remap.rep_igl.size()),
                                                                     remap.rep_igl.data(),
                                                                     remap.fft_isz.data(),
                                                                     remap.phase.data());
                            }
                            else
                            {
                                wfcpw->recip2real_remapped(psi_nk,
                                                           psi_nk_real,
                                                           static_cast<int>(remap.rep_igl.size()),
                                                           remap.rep_igl.data(),
                                                           remap.fft_isz.data(),
                                                           remap.phase.data());
                            }
                        }
                    }
                }
                const Real k_occ = wg_nkb / wk_nk;

                if (PARAM.inp.exx_use_q_tile)
                {
                    for (int q_start = 0; q_start < static_cast<int>(q_points.size()); q_start += q_tile_size)
                    {
                        const int q_count = std::min(q_tile_size, static_cast<int>(q_points.size()) - q_start);
                        for (int m_start = 0; m_start < nbands_psi; m_start += source_tile_size)
                        {
                            const int m_count = std::min(source_tile_size, nbands_psi - m_start);
                            setmem_complex_op()(q_real_tile,
                                                0,
                                                static_cast<std::size_t>(q_tile_size)
                                                    * static_cast<std::size_t>(source_tile_size)
                                                    * static_cast<std::size_t>(wfcpw->nrxx));
                            std::fill(q_weights.begin(), q_weights.end(), Real(0));

                            for (int q_local = 0; q_local < q_count; ++q_local)
                            {
                                const auto* qpoint = q_points[q_start + q_local];
                                const int iq_rep_spin = p_kv->exx_rep_spin_index(*qpoint, ispin);
                                const bool own_qpoint = qpoint->rep_pool == GlobalV::MY_POOL;
                                for (int m_local = 0; m_local < m_count; ++m_local)
                                {
                                    const int mband = m_start + m_local;
                                    double wg_mqb = 0.0;
                                    double wk_mq = 0.0;
                                    if (own_qpoint)
                                    {
                                        wg_mqb = wg(iq_rep_spin, mband);
                                        wk_mq = p_kv->wk[iq_rep_spin];
                                    }
#ifdef __MPI
                                    if (GlobalV::KPAR > 1)
                                    {
                                        MPI_Bcast(&wg_mqb,
                                                  1,
                                                  MPI_DOUBLE,
                                                  p_kv->para_k.get_startpro_pool(qpoint->rep_pool),
                                                  MPI_COMM_WORLD);
                                        MPI_Bcast(&wk_mq,
                                                  1,
                                                  MPI_DOUBLE,
                                                  p_kv->para_k.get_startpro_pool(qpoint->rep_pool),
                                                  MPI_COMM_WORLD);
                                    }
#endif
                                    const std::size_t tile_state = static_cast<std::size_t>(q_local) * source_tile_size
                                                                   + m_local;
                                    if (wg_mqb >= 1e-12)
                                    {
                                        q_weights[tile_state] = static_cast<Real>(wg_mqb / wk_mq * qpoint->weight);
                                        if (own_qpoint)
                                        {
                                            load_q_real(*qpoint,
                                                        ispin,
                                                        mband,
                                                        q_real_tile
                                                            + tile_state * static_cast<std::size_t>(wfcpw->nrxx));
                                        }
#ifdef __MPI
                                        if (GlobalV::KPAR > 1)
                                        {
                                            Parallel_Common::bcast_dev<T, Device>(
                                                q_real_tile + tile_state * static_cast<std::size_t>(wfcpw->nrxx),
                                                wfcpw->nrxx,
                                                KP_WORLD,
                                                qpoint->rep_pool);
                                        }
#endif
                                    }
                                }
                            }

                            if (!own_kpoint || wg_nkb < 1e-12)
                            {
                                continue;
                            }

                            for (int q_local = 0; q_local < q_count; ++q_local)
                            {
                                const auto* qpoint = q_points[q_start + q_local];
                                hamilt::get_exx_potential<Real, Device>(p_kv,
                                                                        wfcpw,
                                                                        rhopw_exx,
                                                                        pot,
                                                                        tpiba,
                                                                        gamma_extrapolation,
                                                                        omega,
                                                                        kpoint,
                                                                        *qpoint,
                                                                        true);
                                hamilt::get_exx_stress_potential<Real, Device>(p_kv,
                                                                               wfcpw,
                                                                               rhopw_exx,
                                                                               pot_stress,
                                                                               tpiba,
                                                                               gamma_extrapolation,
                                                                               omega,
                                                                               kpoint,
                                                                               *qpoint);

                                for (int m_local = 0; m_local < m_count; ++m_local)
                                {
                                    const Real q_weight = q_weights[static_cast<std::size_t>(q_local) * source_tile_size
                                                                    + m_local];
                                    if (std::abs(q_weight) < 1e-12)
                                    {
                                        continue;
                                    }
                                    T* psi_mq_real_tile
                                        = q_real_tile
                                          + (static_cast<std::size_t>(q_local) * source_tile_size + m_local)
                                                * static_cast<std::size_t>(wfcpw->nrxx);

                        // overlap density in real space
                        setmem_complex_op()(density_real, 0.0, rhopw_exx->nrxx);
                        hamilt::cal_density_real_op<T, Device>()(psi_nk_real,
                                                                  psi_mq_real_tile,
                                                                  density_real,
                                                                  omega,
                                                                  wfcpw->nrxx);

                        // density in reciprocal space
                        rhopw_exx->template real_to_recip<T, T, Device>(density_real, density_recip);

                        // really calculate the stress
                        if (use_gpu_q_tile)
                        {
                            const Real scalar = static_cast<Real>(-GlobalC::exx_info.info_global.hybrid_alpha
                                                                  * 0.25 * k_occ * kpoint.weight
                                                                  * k_spin_degeneracy * q_weight);
                            hamilt::exx_stress_accumulate_op<T, Device>()(
                                density_recip,
                                pot,
                                pot_stress,
                                gcar_flat,
                                static_cast<Real>(kpoint.full_kvec_c.x - qpoint->full_kvec_c.x),
                                static_cast<Real>(kpoint.full_kvec_c.y - qpoint->full_kvec_c.y),
                                static_cast<Real>(kpoint.full_kvec_c.z - qpoint->full_kvec_c.z),
                                static_cast<Real>(tpiba),
                                scalar,
                                rhopw_exx->npw,
                                sigma_exx_device);
                            continue;
                        }

                        // for alpha beta
                        for (int alpha = 0; alpha < 3; alpha++)
                        {
                            for (int beta = alpha; beta < 3; beta++)
                            {
                                int delta_ab = (alpha == beta) ? 1 : 0;
                                double sigma_ab_loc = 0.0;
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(static) reduction(+:sigma_ab_loc)
                            #endif
                                for (int ig = 0; ig < rhopw_exx->npw; ig++)
                                {
                                    const ModuleBase::Vector3<double> kqg
                                        = kpoint.full_kvec_c - qpoint->full_kvec_c + rhopw_exx->gcar[ig];
                                    double kqg_alpha = kqg[alpha] * tpiba;
                                    double kqg_beta = kqg[beta] * tpiba;
                                    // equation 10 of 10.1103/PhysRevB.73.125120
                                    double density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
                                    const int idx = ig;
                                    double pot_local = pot[idx];
                                    double pot_stress_local = pot_stress[idx];
                                    sigma_ab_loc += density_recip2 * pot_local
                                                    * (kqg_alpha * kqg_beta * pot_stress_local - delta_ab);

                                }

                                // 0.5 in the following line is caused by 2x in the pot
                                sigma(alpha, beta) -= GlobalC::exx_info.info_global.hybrid_alpha
                                                      * 0.25 * sigma_ab_loc * k_occ * kpoint.weight
                                                      * k_spin_degeneracy * q_weight;
                            }
                        }
                                }
                            }
                        }
                    }
                    continue;
                }

                for (const auto& qpoint: p_kv->exx_full_q_map)
                {
                    if (!qpoint.active)
                    {
                        continue;
                    }
                    const int iq_rep_spin = p_kv->exx_rep_spin_index(qpoint, ispin);
                    if (own_kpoint && wg_nkb >= 1e-12)
                    {
                        hamilt::get_exx_potential<Real, Device>(p_kv,
                                                                wfcpw,
                                                                rhopw,
                                                                pot,
                                                                tpiba,
                                                                gamma_extrapolation,
                                                                omega,
                                                                kpoint,
                                                                qpoint,
                                                                true);
                        hamilt::get_exx_stress_potential<Real, Device>(p_kv,
                                                                       wfcpw,
                                                                       rhopw,
                                                                       pot_stress,
                                                                       tpiba,
                                                                       gamma_extrapolation,
                                                                       omega,
                                                                       kpoint,
                                                                       qpoint);
                    }
                    const bool own_qpoint = qpoint.rep_pool == GlobalV::MY_POOL;
                    for (int mband = 0; mband < d_psi_in->get_nbands(); mband++)
                    {
                        double wg_mqb = 0.0;
                        double wk_mq = 0.0;
                        if (own_qpoint)
                        {
                            wg_mqb = wg(iq_rep_spin, mband);
                            wk_mq = p_kv->wk[iq_rep_spin];
                        }
#ifdef __MPI
                        MPI_Bcast(&wg_mqb, 1, MPI_DOUBLE, p_kv->para_k.get_startpro_pool(qpoint.rep_pool), MPI_COMM_WORLD);
                        MPI_Bcast(&wk_mq, 1, MPI_DOUBLE, p_kv->para_k.get_startpro_pool(qpoint.rep_pool), MPI_COMM_WORLD);
#endif
                        if (wg_mqb < 1e-12 || wg_nkb < 1e-12)
                        {
                            continue;
                        }
                        if (own_qpoint)
                        {
                            load_q_real(qpoint, ispin, mband, psi_mq_real);
                        }
#ifdef __MPI
                        Parallel_Common::bcast_dev<T, Device>(psi_mq_real, wfcpw->nrxx, KP_WORLD, qpoint.rep_pool);
#endif
                        if (!own_kpoint)
                        {
                            continue;
                        }

                        // overlap density in real space
                        setmem_complex_op()(density_real, 0.0, rhopw->nrxx);
                        for (int ig = 0; ig < rhopw->nrxx; ig++)
                        {
                            density_real[ig] = psi_nk_real[ig] * std::conj(psi_mq_real[ig]) * omega_inv;
                        }

                        // density in reciprocal space
                        rhopw->template real_to_recip<T, T, Device>(density_real, density_recip);

                        // really calculate the stress

                        // for alpha beta
                        for (int alpha = 0; alpha < 3; alpha++)
                        {
                            for (int beta = alpha; beta < 3; beta++)
                            {
                                int delta_ab = (alpha == beta) ? 1 : 0;
                                double sigma_ab_loc = 0.0;
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(static) reduction(+:sigma_ab_loc)
                            #endif
                                for (int ig = 0; ig < rhopw->npw; ig++)
                                {
                                    const ModuleBase::Vector3<double> kqg
                                        = kpoint.full_kvec_c - qpoint.full_kvec_c + rhopw->gcar[ig];
                                    double kqg_alpha = kqg[alpha] * tpiba;
                                    double kqg_beta = kqg[beta] * tpiba;
                                    // equation 10 of 10.1103/PhysRevB.73.125120
                                    double density_recip2 = std::real(density_recip[ig] * std::conj(density_recip[ig]));
                                    const int idx = ig;
                                    double pot_local = pot[idx];
                                    double pot_stress_local = pot_stress[idx];
                                    sigma_ab_loc += density_recip2 * pot_local
                                                    * (kqg_alpha * kqg_beta * pot_stress_local - delta_ab);

                                }

                                // 0.5 in the following line is caused by 2x in the pot
                                const double q_occ = wg_mqb / wk_mq;
                                sigma(alpha, beta) -= GlobalC::exx_info.info_global.hybrid_alpha
                                                      * 0.25 * sigma_ab_loc * k_occ * kpoint.weight
                                                      * k_spin_degeneracy * q_occ * qpoint.weight;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int l = 0; l < 3; l++)
    {
        for (int m = l + 1; m < 3; m++)
        {
            sigma(m, l) = sigma(l, m);
        }
    }
    if (use_gpu_q_tile)
    {
        Real sigma_exx_host[6] = {0, 0, 0, 0, 0, 0};
        syncmem_real_d2h_op()(sigma_exx_host, sigma_exx_device, 6);
        int idx = 0;
        for (int alpha = 0; alpha < 3; ++alpha)
        {
            for (int beta = alpha; beta < 3; ++beta)
            {
                sigma(alpha, beta) += sigma_exx_host[idx++];
            }
        }
        for (int l = 0; l < 3; l++)
        {
            for (int m = l + 1; m < 3; m++)
            {
                sigma(m, l) = sigma(l, m);
            }
        }
    }

    // Full-k terms are accumulated only on their representative owning pool,
    // while each pool's FFT/G work remains distributed over ranks.
    Parallel_Reduce::reduce_all(sigma.c, sigma.nr * sigma.nc);


    delmem_complex_op()(psi_nk_real);
    delmem_complex_op()(psi_mq_real);
    delmem_complex_op()(density_real);
    delmem_complex_op()(density_recip);
    delmem_real_op()(pot);
    delmem_real_op()(pot_stress);
    delmem_real_op()(sigma_exx_device);
    delmem_real_op()(gcar_flat);
    delmem_complex_op()(q_real_tile);
    delete rhopw_exx_owned;
}

template class Stress_PW<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_PW<double, base_device::DEVICE_GPU>;
#endif

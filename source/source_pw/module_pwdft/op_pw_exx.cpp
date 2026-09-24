#include "op_pw_exx.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_comm.h" // use KP_WORLD
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/klist.h"
#include "source_hamilt/operator.h"
#include "source_psi/psi.h"
#include "source_pw/module_pwdft/kernels/cal_density_real_op.h"
#include "source_pw/module_pwdft/kernels/exx_batch_op.h"
#include "source_pw/module_pwdft/kernels/exx_cal_energy_op.h"
#include "source_pw/module_pwdft/kernels/mul_potential_op.h"
#include "source_pw/module_pwdft/kernels/vec_mul_cx_op.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <utility>

namespace hamilt
{
template <typename T, typename Device>
std::vector<typename GetTypeReal<T>::type> OperatorEXXPW<T, Device>::fock_div = {};

template <typename T, typename Device>
std::vector<typename GetTypeReal<T>::type> OperatorEXXPW<T, Device>::erfc_div = {};

template <typename T, typename Device>
OperatorEXXPW<T, Device>::OperatorEXXPW(const int* isk_in,
                                        const ModulePW::PW_Basis_K* wfcpw_in,
                                        const ModulePW::PW_Basis* rhopw_in,
                                        K_Vectors *kv_in,
                                        const UnitCell *ucell,
                                        const bool separate_loop_in,
                                        const Real hybrid_alpha_in,
                                        const CoulombParam& coulomb_param_in)
    : isk(isk_in), wfcpw(wfcpw_in), rhopw(rhopw_in), kv(kv_in), ucell(ucell),
      separate_loop(separate_loop_in), hybrid_alpha(hybrid_alpha_in),
      coulomb_param(coulomb_param_in)
{
    if (GlobalV::KPAR != 1 && PARAM.inp.exxace == false)
    {
        // GlobalV::ofs_running << "EXX Calculation does not support k-point parallelism" << std::endl;
        ModuleBase::WARNING_QUIT("OperatorEXXPW", "EXX Calculation does not support k-point parallelism when exxace is set to false");
    }
    gamma_extrapolation = PARAM.inp.exx_gamma_extrapolation;
    bool is_mp = kv_in->get_is_mp();
#ifdef __MPI
    Parallel_Common::bcast_bool(is_mp);
#endif
    if (!is_mp)
    {
        gamma_extrapolation = false;
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
    resmem_real_op()(pot, rhopw->npw);

    tpiba = ucell->tpiba;
    Real tpiba2 = tpiba * tpiba;

    // initialize rhopw_dev
    double ecut_exx = PARAM.inp.ecutexx;
    if (ecut_exx == 0.0)
    {
        ecut_exx = PARAM.inp.ecutrho;
    }

    rhopw_dev = new ModulePW::PW_Basis(wfcpw->get_device(), rhopw->get_precision());
    rhopw_dev->fft_bundle.setfft(wfcpw->get_device(), rhopw->get_precision());
#ifdef __MPI
    rhopw_dev->initmpi(rhopw->poolnproc, rhopw->poolrank, rhopw->pool_world);
#endif
    // here we can actually use different ecut to init the grids
    rhopw_dev->initgrids(rhopw->lat0, rhopw->latvec, ecut_exx);
    rhopw_dev->initgrids(rhopw->lat0, rhopw->latvec, rhopw->nx, rhopw->ny, rhopw->nz);
    rhopw_dev->initparameters(rhopw->gamma_only, ecut_exx, rhopw->distribution_type, rhopw->xprime);
    rhopw_dev->setuptransform();
    rhopw_dev->collect_local_pw();

    auto param_fock = this->coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock];
    for (auto param: param_fock)
    {
        fock_div.push_back(exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type::Fock,
                                          0.0,
                                          kv,
                                          wfcpw,
                                          rhopw_dev,
                                          tpiba,
                                          gamma_extrapolation,
                                          ucell->omega));
    }
    auto param_erfc = this->coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Erfc];
    for (auto param: param_erfc)
    {
        erfc_div.push_back(exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type::Erfc,
                                          std::stod(param["omega"]),
                                          kv,
                                          wfcpw,
                                          rhopw_dev,
                                          tpiba,
                                          gamma_extrapolation,
                                          ucell->omega));
    }

}   // end of constructor

template <typename T, typename Device>
OperatorEXXPW<T, Device>::~OperatorEXXPW()
{
    // use delete_memory_op to delete the allocated pws
    delmem_complex_op()(psi_nk_real);
    delmem_complex_op()(psi_mq_real);
    if (psi_nk_real_cache != nullptr)
    {
        delmem_complex_op()(psi_nk_real_cache);
    }
    if (dens_box_batch != nullptr)
    {
        delmem_complex_op()(dens_box_batch);
        delmem_complex_op()(dens_pw_batch);
        exx_batch_fft_plan_destroy<T, Device>(&exx_fft_plan);
        exx_batch_fft_plan_destroy<T, Device>(&exx_fft_plan1);
    }
    if (sg_map_rho != nullptr)
    {
        delmem_int_op()(sg_map_rho);
    }
    if (sg_map_wfc != nullptr)
    {
        delmem_int_op()(sg_map_wfc);
    }
    if (full_map_rho != nullptr)
    {
        delmem_int_op()(full_map_rho);
    }
    if (full_map_wfc != nullptr)
    {
        delmem_int_op()(full_map_wfc);
    }
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
    // std::cout << cal_exx_energy_ace(&psi) << " EXX energy" << std::endl;
    // MPI_Abort(MPI_COMM_WORLD, 0);
    // return;

    if (is_first_node)
    {
        setmem_complex_op()(tmhpsi, 0, nbasis*nbands/npol);
    }

    if (PARAM.inp.exxace && this->separate_loop)
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
    ModuleBase::timer::start("OperatorEXXPW", "act_op");

    setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
    setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_recip, 0, rhopw_dev->npw);
    setmem_complex_op()(psi_nk_real, 0, wfcpw->nrxx);
    setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);

    auto q_points = get_q_points(this->ik);
    int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    int nk = wfcpw->nks / nk_fac;
    const Real nqs = q_points.size();

    maybe_setup_exx_grid();
    // psi_nk in real space for all bands once per ik, reused over (iq, m)
    cache_psi_nk_real(nbands, nbasis, tmpsi_in, this->ik);

    for (int iq: q_points)
    {
        get_exx_potential<Real, Device>(kv, wfcpw, rhopw_dev, pot, tpiba, gamma_extrapolation, ucell->omega, this->ik, iq % nk, false, this->coulomb_param);
        for (int m_iband = 0; m_iband < psi.get_nbands(); m_iband++)
        {
            // occupation of the source state (m, iq), not of the target k-point
            double wg_mqb_real = (*wg)(iq, m_iband);
            if (wg_mqb_real < 1e-12)
            {
                continue;
            }

            wfc_to_real_exx(get_pw(m_iband, iq), iq, nbasis);

            // full accumulation weight, hybrid_alpha included
            const Real factor = this->hybrid_alpha * wg_mqb_real / kv->wk[iq] / nqs;
            apply_fock_all_bands(nbands, nbasis, iq, factor, tmhpsi);

        } // end of m_iband

    } // end of iq

    ModuleBase::timer::end("OperatorEXXPW", "act_op");

}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::act_op_kpar(const int nbands,
                                   const int nbasis,
                                   const int npol,
                                   const T *tmpsi_in,
                                   T *tmhpsi,
                                   const int ngk_ik,
                                   const bool is_first_node) const
{
    ModuleBase::timer::start("OperatorEXXPW", "act_op_kpar");

    setmem_complex_op()(h_psi_recip, 0, wfcpw->npwk_max);
    setmem_complex_op()(h_psi_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_real, 0, rhopw_dev->nrxx);
    setmem_complex_op()(density_recip, 0, rhopw_dev->npw);
    // setmem_complex_op()(psi_all_real, 0, wfcpw->nrxx * GlobalV::NBANDS);
    // std::map<std::pair<int, int>, bool> has_real;
    setmem_complex_op()(psi_nk_real, 0, wfcpw->nrxx);
    setmem_complex_op()(psi_mq_real, 0, wfcpw->nrxx);
    int nqs = kv->get_nkstot_nospin();
    int nspin_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    int ispin = this->ik < (wfcpw->nks / nspin_fac) ? 0 : 1;

    maybe_setup_exx_grid();
    // psi_nk in real space for all bands once per ik, reused over (iq, m);
    // the MPI communication order below is unchanged
    cache_psi_nk_real(nbands, nbasis, tmpsi_in, this->ik);

    // ik fixed here, select band n
    for (int iq = 0; iq < nqs; iq++)
    {
        // for \psi_nk, get the pw of iq and band m
        get_exx_potential<Real,  Device>(kv, wfcpw, rhopw_dev, pot, tpiba, gamma_extrapolation, ucell->omega, this->ik, iq, false, this->coulomb_param);

        // decide which pool does the iq belong to
        int iq_pool = kv->para_k.whichpool[iq];
        int iq_loc  = iq - kv->para_k.startk_pool[iq_pool];
        int iq_loc_spin = iq_loc;
        if (ispin == 1)
        {
            iq_loc_spin += wfcpw->nks / nspin_fac;
        }

        // occupation row and k weight of the source k-point, fetched from the
        // pool that owns it in a single broadcast
        const int nb = psi.get_nbands();
        std::vector<double> occ_q(nb + 1);
        if (iq_pool == GlobalV::MY_POOL)
        {
            for (int m = 0; m < nb; m++)
            {
                occ_q[m] = (*wg)(iq_loc_spin, m);
            }
            occ_q[nb] = kv->wk[iq_loc_spin];
        }
#ifdef __MPI
        Parallel_Common::bcast_dev<double, base_device::DEVICE_CPU>(
            occ_q.data(), nb + 1, MPI_COMM_WORLD, kv->para_k.get_startpro_pool(iq_pool));
#endif
        const Real wk_q = occ_q[nb];

        for (int m_iband = 0; m_iband < nb; m_iband++)
        {
            const double wg_mqb = occ_q[m_iband];
            if (wg_mqb < 1e-12)
                continue;

            if (iq_pool == GlobalV::MY_POOL)
            {
                wfc_to_real_exx(get_pw(m_iband, iq_loc_spin), iq_loc, nbasis);
            }
#ifdef __MPI
            Parallel_Common::bcast_dev<T, Device>(psi_mq_real, exx_grid_size(), KP_WORLD, iq_pool);
#endif

            // k weight of the source k-point (identical to wk[this->ik] on the
            // uniform k-grids without symmetry reduction this scheme assumes)
            const Real factor = this->hybrid_alpha * wg_mqb / wk_q / nqs;
            apply_fock_all_bands(nbands, nbasis, iq, factor, tmhpsi);

        } // end of m_iband

    } // end of iq

    ModuleBase::timer::end("OperatorEXXPW", "act_op_kpar");

}


template <typename T, typename Device>
void OperatorEXXPW<T, Device>::maybe_setup_exx_grid() const
{
    if (exx_sg_init)
    {
        return;
    }
#if !defined(__CUDA)
    // the small-grid/batched path on GPU is implemented for CUDA only; on
    // ROCm builds the GPU instantiation keeps the full-grid path
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        exx_sg_init = true;
        return;
    }
#endif
    setup_exx_small_grid();
    if (std::is_same<Device, base_device::DEVICE_CPU>::value && !exx_sg_ok)
    {
        // the GPU bases own these maps as ig2ixyz*; the CPU batched path
        // needs host-built ones for the full grid
        setup_full_grid_maps();
    }
}

template <typename T, typename Device>
bool OperatorEXXPW<T, Device>::batch_active() const
{
#if !defined(__CUDA)
    if (!std::is_same<Device, base_device::DEVICE_CPU>::value)
    {
        return false; // ROCm GPU
    }
#endif
    // batched FFTs need the whole box local to this rank; on CPU the small
    // grid additionally cuts the FFT work, the full-grid batch mainly keeps
    // the code path unified
    return exx_sg_ok
           || (wfcpw->nrxx == rhopw_dev->nrxx && wfcpw->nxyz == rhopw_dev->nxyz
               && wfcpw->nrxx == wfcpw->nxyz);
}

template <typename T, typename Device>
const int* OperatorEXXPW<T, Device>::active_map_rho() const
{
    if (exx_sg_ok)
    {
        return sg_map_rho;
    }
#if defined(__CUDA)
    if (std::is_same<Device, base_device::DEVICE_GPU>::value)
    {
        return rhopw_dev->ig2ixyz_gpu;
    }
#endif
    return full_map_rho;
}

template <typename T, typename Device>
const int* OperatorEXXPW<T, Device>::active_map_wfc(const int ik) const
{
    if (exx_sg_ok)
    {
        return sg_map_wfc + ik * wfcpw->npwk_max;
    }
#if defined(__CUDA)
    if (std::is_same<Device, base_device::DEVICE_GPU>::value)
    {
        return wfcpw->ig2ixyz_k + ik * wfcpw->npwk_max;
    }
#endif
    return full_map_wfc + ik * wfcpw->npwk_max;
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::cache_psi_nk_real(const int nbands, const int nbasis, const T* psi_in, const int ik) const
{
    if (psi_nk_cache_size < nbands)
    {
        // resize_memory_op frees the old block itself
        resmem_complex_op()(psi_nk_real_cache, nbands * exx_grid_size());
        psi_nk_cache_size = nbands;
    }
    if (batch_active())
    {
        // scatter the PW coefficients into the active box and do one batched
        // backward FFT (unnormalized, like recip_to_real)
        ensure_exx_batch(nbands);
        setmem_complex_op()(psi_nk_real_cache, 0, nbands * exx_grid_size());
        exx_batch_scatter_wfc<T, Device>(nbands,
                                         wfcpw->npwk[ik],
                                         exx_grid_size(),
                                         active_map_wfc(ik),
                                         psi_in,
                                         nbasis,
                                         psi_nk_real_cache);
        exx_batch_fft_exec<T, Device>(exx_fft_plan, psi_nk_real_cache, false);
        return;
    }
    for (int n = 0; n < nbands; n++)
    {
        wfcpw->recip_to_real(ctx, psi_in + n * nbasis, psi_nk_real_cache + n * wfcpw->nrxx, ik);
    }
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::wfc_to_real_exx(const T* psi_m, const int iq, const int psi_stride) const
{
    if (batch_active())
    {
        // scatter + backward FFT on the active grid
        setmem_complex_op()(psi_mq_real, 0, exx_grid_size());
        exx_batch_scatter_wfc<T, Device>(1,
                                         wfcpw->npwk[iq],
                                         exx_grid_size(),
                                         active_map_wfc(iq),
                                         psi_m,
                                         psi_stride,
                                         psi_mq_real);
        exx_batch_fft_exec<T, Device>(exx_fft_plan1, psi_mq_real, false);
        return;
    }
    wfcpw->recip_to_real(ctx, psi_m, psi_mq_real, iq);
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::apply_fock_all_bands(const int nbands,
                                                    const int nbasis,
                                                    const int iq,
                                                    const Real factor,
                                                    T* tmhpsi) const
{
    if (batch_active())
    {
        ensure_exx_batch(nbands);
        apply_exx_nbatched(nbands, nbasis, psi_mq_real, factor, tmhpsi);
        return;
    }
    for (int n = 0; n < nbands; n++)
    {
        const T* nk_real = psi_nk_real_cache + n * wfcpw->nrxx;

        // direct multiplication in real space, \psi_nk(r) * \psi_mq(r)
        cal_density_recip(nk_real, psi_mq_real, ucell->omega);

        // multiply the density with the potential in recip space
        multiply_potential(density_recip, this->ik, iq);

        // bring the potential back to real space
        rho_recip2real(density_recip, density_real);

        vec_mul_vec_complex_op<T, Device>()(density_real, psi_mq_real, density_real, wfcpw->nrxx);

        wfcpw->real_to_recip(ctx, density_real, tmhpsi + n * nbasis, this->ik, true, factor);
    }
}


template <typename T, typename Device>
void OperatorEXXPW<T, Device>::prepare_pair_densities(const int nbands) const
{
    if (batch_active())
    {
        // pair densities of all bands in one batched round, rhopw_dev G-space
        ensure_exx_batch(nbands);
        calc_density_pw_nbatched(nbands, psi_mq_real);
    }
}

template <typename T, typename Device>
const T* OperatorEXXPW<T, Device>::pair_density(const int n) const
{
    if (batch_active())
    {
        return dens_pw_batch + n * rhopw_dev->npw;
    }
    cal_density_recip(psi_nk_real_cache + n * wfcpw->nrxx, psi_mq_real, ucell->omega);
    return density_recip;
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::setup_exx_small_grid() const
{
    exx_sg_init = true;
    exx_sg_ok = false;

    const bool user_set = PARAM.inp.ecutexx > 0.0;
    double ecut_exx = user_set ? PARAM.inp.ecutexx : PARAM.inp.ecutrho;

    // FFT box dims for ecut_exx, obtained the same way as rhopw_dev
    ModulePW::PW_Basis gridt(wfcpw->get_device(), rhopw->get_precision());
    gridt.fft_bundle.setfft(wfcpw->get_device(), rhopw->get_precision());
#ifdef __MPI
    gridt.initmpi(rhopw->poolnproc, rhopw->poolrank, rhopw->pool_world);
#endif
    gridt.initgrids(rhopw->lat0, rhopw->latvec, ecut_exx);
    if (gridt.nx == wfcpw->nx && gridt.ny == wfcpw->ny && gridt.nz == wfcpw->nz)
    {
        if (user_set && GlobalV::MY_RANK == 0)
        {
            ModuleBase::WARNING("OperatorEXXPW",
                                "ecutexx gives no smaller FFT grid than ecutrho; EXX stays on the full grid");
        }
        return; // no gain over the current grid
    }
    if (rhopw->poolnproc != 1)
    {
        if (user_set && GlobalV::MY_RANK == 0)
        {
            ModuleBase::WARNING("OperatorEXXPW",
                                "ecutexx is set but the FFT box is distributed over the plane-wave pool; "
                                "the small-grid EXX path needs a local box, so EXX stays on the full grid");
        }
        return; // batched path assumes the whole box is local
    }

    // QE-style guard (gcutmt >= gkcut): every |k+G|^2 must fit the ecut_exx
    // sphere, otherwise the wavefunctions cannot live on the small grid.
    // Note: gk2 entries beyond npwk[ik] are uninitialized, so loop per k.
    const Real tpiba2 = tpiba * tpiba;
    for (int ik = 0; ik < wfcpw->nks; ik++)
    {
        for (int ig = 0; ig < wfcpw->npwk[ik]; ig++)
        {
            if (wfcpw->gk2[ik * wfcpw->npwk_max + ig] * tpiba2 > ecut_exx)
            {
                if (user_set && GlobalV::MY_RANK == 0)
                {
                    ModuleBase::WARNING("OperatorEXXPW",
                                        "ecutexx is smaller than the wavefunction cutoff |k+G|^2 of some "
                                        "k-point; the wavefunctions would not fit the small grid, so EXX "
                                        "stays on the full grid. Raise ecutexx above ~ecutwfc to enable it");
                }
                return;
            }
        }
    }

    sg_nx = gridt.nx;
    sg_ny = gridt.ny;
    sg_nz = gridt.nz;
    sg_nxyz = sg_nx * sg_ny * sg_nz;

    // Remap a box index of the (nx, ny, nz) grid to the small grid, going
    // through the Miller indices (both boxes are centered at G=0).
    auto remap = [](const int idx, const int nx, const int ny, const int nz,
                    const int tnx, const int tny, const int tnz) {
        const int iz = idx % nz;
        const int iy = (idx / nz) % ny;
        const int ix = idx / (ny * nz);
        const int hx = ix >= nx / 2 + 1 ? ix - nx : ix;
        const int hy = iy >= ny / 2 + 1 ? iy - ny : iy;
        const int hz = iz >= nz / 2 + 1 ? iz - nz : iz;
        const int tx = hx < 0 ? hx + tnx : hx;
        const int ty = hy < 0 ? hy + tny : hy;
        const int tz = hz < 0 ? hz + tnz : hz;
        return tz + ty * tnz + tx * tny * tnz;
    };

    // rho map: rhopw_dev G-vectors (box = wfcpw/rhopw grid) -> small box
    const int npw_rho = rhopw_dev->npw;
    std::vector<int> host_map(npw_rho);
    std::vector<int> host_ig2ixyz(npw_rho);
    if (rhopw_dev->ig2ixyz_gpu != nullptr)
    {
        syncmem_int_d2h_op()(host_ig2ixyz.data(), rhopw_dev->ig2ixyz_gpu, npw_rho);
    }
    else
    {
        // host basis without ig2ixyz: rebuild the big-box index from the
        // Miller indices (gdirect holds integer values as doubles)
        for (int ig = 0; ig < npw_rho; ig++)
        {
            const int hx = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].x));
            const int hy = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].y));
            const int hz = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].z));
            const int ix = hx < 0 ? hx + wfcpw->nx : hx;
            const int iy = hy < 0 ? hy + wfcpw->ny : hy;
            const int iz = hz < 0 ? hz + wfcpw->nz : hz;
            host_ig2ixyz[ig] = iz + iy * wfcpw->nz + ix * wfcpw->ny * wfcpw->nz;
        }
    }
    for (int ig = 0; ig < npw_rho; ig++)
    {
        host_map[ig] = remap(host_ig2ixyz[ig], wfcpw->nx, wfcpw->ny, wfcpw->nz, sg_nx, sg_ny, sg_nz);
    }
    resmem_int_op()(sg_map_rho, npw_rho);
    syncmem_int_h2d_op()(sg_map_rho, host_map.data(), npw_rho);

    // wfc map: (G+k) components of every k-point -> small box
    const int n_tot = wfcpw->npwk_max * wfcpw->nks;
    host_map.resize(n_tot);
    if (!wfcpw->ig2ixyz_k_cpu.empty())
    {
        for (int i = 0; i < n_tot; i++)
        {
            host_map[i] = remap(wfcpw->ig2ixyz_k_cpu[i], wfcpw->nx, wfcpw->ny, wfcpw->nz, sg_nx, sg_ny, sg_nz);
        }
    }
    else
    {
        // host basis without ig2ixyz_k_cpu: same construction as
        // PW_Basis_K::get_ig2ixyz_k, from the stick decomposition (the small
        // grid requires a local box, so all sticks are on this rank)
        for (int ik = 0; ik < wfcpw->nks; ik++)
        {
            for (int igl = 0; igl < wfcpw->npwk[ik]; igl++)
            {
                const int isz = wfcpw->igl2isz_k[igl + ik * wfcpw->npwk_max];
                const int iz = isz % wfcpw->nz;
                const int is = isz / wfcpw->nz;
                const int ixy = wfcpw->is2fftixy[is];
                const int iy = ixy % wfcpw->ny;
                const int ix = ixy / wfcpw->ny;
                host_map[igl + ik * wfcpw->npwk_max] =
                    remap(iz + iy * wfcpw->nz + ix * wfcpw->ny * wfcpw->nz,
                          wfcpw->nx, wfcpw->ny, wfcpw->nz, sg_nx, sg_ny, sg_nz);
            }
        }
    }
    resmem_int_op()(sg_map_wfc, n_tot);
    syncmem_int_h2d_op()(sg_map_wfc, host_map.data(), n_tot);

    exx_sg_ok = true;
    if (GlobalV::MY_RANK == 0)
    {
        std::cout << " EXX small grid enabled: FFT (" << sg_nx << "," << sg_ny << "," << sg_nz << ") instead of ("
                  << wfcpw->nx << "," << wfcpw->ny << "," << wfcpw->nz << ") for ecut_exx = " << ecut_exx << " Ry"
                  << std::endl;
    }
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::setup_full_grid_maps() const
{
    // host-built full-box maps for the CPU batched path when the small grid
    // is not active; same sources as in setup_exx_small_grid
    const int npw_rho = rhopw_dev->npw;
    std::vector<int> host_map(npw_rho);
    for (int ig = 0; ig < npw_rho; ig++)
    {
        const int hx = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].x));
        const int hy = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].y));
        const int hz = static_cast<int>(std::lround(rhopw_dev->gdirect[ig].z));
        const int ix = hx < 0 ? hx + wfcpw->nx : hx;
        const int iy = hy < 0 ? hy + wfcpw->ny : hy;
        const int iz = hz < 0 ? hz + wfcpw->nz : hz;
        host_map[ig] = iz + iy * wfcpw->nz + ix * wfcpw->ny * wfcpw->nz;
    }
    resmem_int_op()(full_map_rho, npw_rho);
    syncmem_int_h2d_op()(full_map_rho, host_map.data(), npw_rho);

    const int n_tot = wfcpw->npwk_max * wfcpw->nks;
    host_map.resize(n_tot);
    for (int ik = 0; ik < wfcpw->nks; ik++)
    {
        for (int igl = 0; igl < wfcpw->npwk[ik]; igl++)
        {
            const int isz = wfcpw->igl2isz_k[igl + ik * wfcpw->npwk_max];
            const int iz = isz % wfcpw->nz;
            const int is = isz / wfcpw->nz;
            const int ixy = wfcpw->is2fftixy[is];
            const int iy = ixy % wfcpw->ny;
            const int ix = ixy / wfcpw->ny;
            host_map[igl + ik * wfcpw->npwk_max] = iz + iy * wfcpw->nz + ix * wfcpw->ny * wfcpw->nz;
        }
    }
    resmem_int_op()(full_map_wfc, n_tot);
    syncmem_int_h2d_op()(full_map_wfc, host_map.data(), n_tot);
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::ensure_exx_batch(const int nbands) const
{
    const int nx = exx_sg_ok ? sg_nx : wfcpw->nx;
    const int ny = exx_sg_ok ? sg_ny : wfcpw->ny;
    const int nz = exx_sg_ok ? sg_nz : wfcpw->nz;
    const int nxyz = exx_sg_ok ? sg_nxyz : wfcpw->nxyz;
    if (exx_batch_alloc == nbands && exx_fft_plan_nx == nx && exx_fft_plan_ny == ny && exx_fft_plan_nz == nz)
    {
        return;
    }
    if (dens_box_batch != nullptr)
    {
        // resize_memory_op frees the old block itself; the plans need an
        // explicit destroy (they are plain library handles)
        exx_batch_fft_plan_destroy<T, Device>(&exx_fft_plan);
        exx_batch_fft_plan_destroy<T, Device>(&exx_fft_plan1);
    }
    resmem_complex_op()(dens_box_batch, nbands * nxyz);
    resmem_complex_op()(dens_pw_batch, nbands * rhopw_dev->npw);
    exx_batch_fft_plan_create<T, Device>(&exx_fft_plan, nx, ny, nz, nbands);
    exx_batch_fft_plan_create<T, Device>(&exx_fft_plan1, nx, ny, nz, 1);
    exx_batch_alloc = nbands;
    exx_fft_plan_nx = nx;
    exx_fft_plan_ny = ny;
    exx_fft_plan_nz = nz;
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::apply_exx_nbatched(const int nbands,
                                                  const int nbasis,
                                                  const T* psi_mq_real,
                                                  const Real factor,
                                                  T* tmhpsi) const
{
    ensure_exx_batch(nbands);
    // active FFT grid: the small ecut_exx grid when usable, else the wfcpw grid
    const int nxyz = exx_sg_ok ? sg_nxyz : wfcpw->nxyz;
    const int nrxx = nxyz; // batched path requires the whole box to be local
    const int npw_rho = rhopw_dev->npw;
    const int npwk = wfcpw->npwk[this->ik];
    const int* map_rho = active_map_rho();
    const int* map_wfc = active_map_wfc(this->ik);

    // 1. reciprocal-space density of all bands, ends up in dens_pw_batch
    calc_density_pw_nbatched(nbands, psi_mq_real);
    // 2. multiply by the Coulomb potential in recip space
    exx_batch_mul_pot<T, Real, Device>(nbands, npw_rho, pot, dens_pw_batch);
    // 3. scatter back to the (pre-zeroed) box
    setmem_complex_op()(dens_box_batch, 0, nbands * nxyz);
    exx_batch_scatter_pw<T, Device>(nbands, npw_rho, nxyz, map_rho, dens_pw_batch, dens_box_batch);
    // 4. batched backward FFT (unnormalized, matching recip2real)
    exx_batch_fft_exec<T, Device>(exx_fft_plan, dens_box_batch, false);
    // 5. multiply by psi_mq(r) in real space
    exx_batch_mul_real<T, Device>(nbands, nrxx, psi_mq_real, dens_box_batch);
    // 6. batched forward FFT
    exx_batch_fft_exec<T, Device>(exx_fft_plan, dens_box_batch, true);
    // 7. gather and accumulate into hpsi with the EXX weight
    exx_batch_gather_accum<T, Real, Device>(nbands, npwk, nxyz, map_wfc, dens_box_batch, factor, tmhpsi, nbasis);
}

template <typename T, typename Device>
void OperatorEXXPW<T, Device>::calc_density_pw_nbatched(const int nbands, const T* psi_mq_real) const
{
    // active FFT grid: the small ecut_exx grid when usable, else the wfcpw grid
    const int nxyz = exx_sg_ok ? sg_nxyz : wfcpw->nxyz;
    const int nrxx = nxyz; // batched path requires the whole box to be local
    const int npw_rho = rhopw_dev->npw;
    const int* map_rho = active_map_rho();

    // 1. density_real(r) = psi_nk(r) * conj(psi_mq(r)) / omega for all bands at once
    exx_batch_density_real<T, Device>(nbands, nrxx, psi_nk_real_cache, psi_mq_real, ucell->omega, dens_box_batch);
    // 2. batched forward FFT to the box
    exx_batch_fft_exec<T, Device>(exx_fft_plan, dens_box_batch, true);
    // 3. gather the PW components (with the 1/nxyz of real_to_recip)
    exx_batch_gather_pw<T, Device>(nbands, npw_rho, nxyz, map_rho, dens_box_batch, dens_pw_batch);
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
    ModuleBase::timer::start("OperatorEXXPW", "multiply_potential");
    int npw = rhopw_dev->npw;
    int nks = wfcpw->nks;
    int nk_fac = PARAM.inp.nspin == 2 ? 2 : 1;
    int nk = nks / nk_fac;

    mul_potential_op<T, Device>()(pot, density_recip, npw, nks, ik, iq);

    ModuleBase::timer::end("OperatorEXXPW", "multiply_potential");
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
    if (PARAM.inp.exxace && this->separate_loop)
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
    const int nb = psi.get_nbands();
    const int nbasis = psi_.get_nbasis();
    const int npw = rhopw_dev->npw;
    maybe_setup_exx_grid();

    double Eexx_ik_real = 0.0;
    for (int ik = 0; ik < wfcpw->nks; ik++)
    {
        // NOTE: psi_nk intentionally comes from the member psi (state of the
        // previous set_psi), while psi_mq comes from the argument psi_ —
        // iter_finish computes dexx as the difference of two cal_exx_energy
        // calls around set_psi, so the asymmetry is load-bearing for the EXX
        // outer-loop convergence check.
        psi.fix_kb(ik, 0);
        cache_psi_nk_real(nb, nbasis, psi.get_pointer(), ik);

        // q points of the same spin channel as ik
        std::vector<int> q_points_ik;
        if (PARAM.inp.nspin == 1)
        {
            for (int iq = 0; iq < wfcpw->nks; iq++)
            {
                q_points_ik.push_back(iq);
            }
        }
        else if (PARAM.inp.nspin == 2)
        {
            const int nk = wfcpw->nks / nk_fac;
            const int k_spin = ik / nk;
            for (int iq = k_spin * nk; iq < (k_spin + 1) * nk; iq++)
            {
                q_points_ik.push_back(iq);
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("OperatorEXXPW", "nspin == 4 not supported");
        }
        const double nqs = q_points_ik.size();
        const int nk = wfcpw->nks / nk_fac;

        for (int iq: q_points_ik)
        {
            get_exx_potential<Real, Device>(kv, wfcpw, rhopw_dev, pot, tpiba, gamma_extrapolation, ucell->omega, ik, iq % nk, false, this->coulomb_param);
            for (int m_iband = 0; m_iband < nb; m_iband++)
            {
                const double wg_iqb_real = (*wg)(iq, m_iband);
                if (wg_iqb_real < 1e-12)
                {
                    continue;
                }
                psi_.fix_kb(iq, m_iband);
                wfc_to_real_exx(psi_.get_pointer(), iq, nbasis);

                // pair densities of all bands on the active grid; with the
                // batched path this is one round for all bands, and the
                // energy kernel below works purely in the G-space
                prepare_pair_densities(nb);
                for (int n_iband = 0; n_iband < nb; n_iband++)
                {
                    const double wg_ikb_real = (*wg)(ik, n_iband);
                    if (wg_ikb_real < 1e-12)
                    {
                        continue;
                    }
                    Eexx_ik_real += exx_cal_energy_op<T, Device>()(pair_density(n_iband),
                                                                   pot,
                                                                   wg_iqb_real / nqs * wg_ikb_real / kv->wk[ik],
                                                                   npw);
                }
            } // m_iband
        } // iq

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

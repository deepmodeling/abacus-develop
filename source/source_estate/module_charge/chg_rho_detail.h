#ifndef CHG_RHO_DETAIL_H
#define CHG_RHO_DETAIL_H

// Internal helpers for charge density mixing (mix_rho_recip/mix_rho_real).
// Not part of the public module_charge API: only chg_mix.cpp
// and the charge mixing unit test are expected to include this header.

#include <functional>
#include <complex>

#include "charge.h"
#include "chg_mix_cfg.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/module_mixing/plain_mixing.h"
#include "source_base/tool_quit.h"

namespace module_charge
{
namespace detail
{

/**
 * @brief Create a two-beta mixing functor: mix the first nunit elements with
 *        mixing_beta and the rest (nunit..total) with mixing_beta_mag.
 *        Used for magnetic cases (nspin==2/4) where the charge channel and
 *        the magnetism channels use different betas.
 * @tparam T element type, double (real space) or std::complex<double> (reciprocal)
 * @param total total number of elements
 * @param nunit number of elements in the charge channel
 * @param mixing_beta beta for the charge channel
 * @param mixing_beta_mag beta for the magnetism channel
 * @return mixing functor
 */
template <typename T>
std::function<void(T*, const T*, const T*)> make_twobeta_mix(
    const int total, const int nunit,
    const double mixing_beta, const double mixing_beta_mag)
{
    return [total, nunit, mixing_beta, mixing_beta_mag](T* out, const T* in, const T* sres)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
        for (int i = 0; i < nunit; ++i)
        {
            out[i] = in[i] + mixing_beta * sres[i];
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
        for (int i = nunit; i < total; ++i)
        {
            out[i] = in[i] + mixing_beta_mag * sres[i];
        }
    };
}

/**
 * @brief Pack charge and magnetism into interleaved layout:
 *        out[0..n]   = d0 + d1  (charge channel)
 *        out[n..2n]  = d0 - d1  (magnetism channel)
 * @tparam T double (real space) or std::complex<double> (reciprocal)
 * @param out output buffer, size >= 2*n
 * @param d0 first component (e.g. chr->rho[0] or chr->rhog[0])
 * @param d1 second component
 * @param n number of elements per component
 */
template <typename T>
void pack_rho_mag(T* out, const T* d0, const T* d1, const int n)
{
    if (n < 0)
    {
        ModuleBase::WARNING_QUIT("pack_rho_mag", "n must be >= 0");
    }
    // A rank may own zero real-space grid points (n == 0) when the grid is
    // decomposed across more processes than it has slabs. Its buffers are
    // legitimately null in that case and the loop below performs no access;
    // only a null buffer with a positive n is a genuine bug.
    if (n > 0 && (out == nullptr || d0 == nullptr || d1 == nullptr))
    {
        ModuleBase::WARNING_QUIT("pack_rho_mag", "pointer is null");
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
    for (int i = 0; i < n; ++i)
    {
        out[i] = d0[i] + d1[i];
        out[i + n] = d0[i] - d1[i];
    }
}

/**
 * @brief Unpack interleaved layout back to charge and magnetism components:
 *        d0[i] = 0.5 * (in[i] + in[i+n])
 *        d1[i] = 0.5 * (in[i] - in[i+n])
 * @tparam T double (real space) or std::complex<double> (reciprocal)
 * @param d0 output first component (e.g. chr->rho[0] or chr->rhog[0])
 * @param d1 output second component
 * @param in input buffer, size >= 2*n
 * @param n number of elements per component
 */
template <typename T>
void unpack_rho_mag(T* d0, T* d1, const T* in, const int n)
{
    if (n < 0)
    {
        ModuleBase::WARNING_QUIT("unpack_rho_mag", "n must be >= 0");
    }
    // A rank may own zero real-space grid points (n == 0); see pack_rho_mag.
    if (n > 0 && (d0 == nullptr || d1 == nullptr || in == nullptr))
    {
        ModuleBase::WARNING_QUIT("unpack_rho_mag", "pointer is null");
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
    for (int i = 0; i < n; ++i)
    {
        d0[i] = 0.5 * (in[i] + in[i + n]);
        d1[i] = 0.5 * (in[i] - in[i + n]);
    }
}

} // namespace detail
} // namespace module_charge

#endif // CHG_RHO_DETAIL_H

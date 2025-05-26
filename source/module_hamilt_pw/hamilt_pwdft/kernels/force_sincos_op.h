#ifndef FORCE_SINCOS_OP_H
#define FORCE_SINCOS_OP_H

#include "module_base/module_device/types.h"
#include <complex>

namespace hamilt
{

template <typename FPTYPE, typename Device>
struct cal_force_loc_sincos_op
{
    /// @brief Calculate local pseudopotential forces (sincos loop only)
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nat - number of atoms
    /// @param npw - number of plane waves
    /// @param gcar - G-vector Cartesian coordinates [npw * 3]
    /// @param tau - atomic positions [nat * 3]
    /// @param iat2it - atom to type mapping [nat]
    /// @param vloc_factors - precomputed vloc factors [npw]
    /// @param aux - charge density in G-space [npw]
    /// @param scale_factor - tpiba * omega
    ///
    /// Output Parameters
    /// @param force - output forces [nat * 3]
    void operator()(const Device* ctx,
                    const int& nat,
                    const int& npw,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    const FPTYPE& scale_factor,
                    FPTYPE* force);
};

template <typename FPTYPE, typename Device>
struct cal_force_ew_sincos_op
{
    /// @brief Calculate Ewald forces (sincos loop only)
    ///
    /// Input Parameters
    /// @param ctx - which device this function runs on
    /// @param nat - number of atoms
    /// @param npw - number of plane waves
    /// @param ig_gge0 - index of G=0 vector (-1 if not present)
    /// @param gcar - G-vector Cartesian coordinates [npw * 3]
    /// @param tau - atomic positions [nat * 3]
    /// @param iat2it - atom to type mapping [nat]
    /// @param it_facts - precomputed it_fact for each atom [nat]
    /// @param aux - structure factor related array [npw]
    ///
    /// Output Parameters
    /// @param force - output forces [nat * 3]
    void operator()(const Device* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

template <typename FPTYPE>
struct cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nat,
                    const int& npw,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    const FPTYPE& scale_factor,
                    FPTYPE* force);
};

template <typename FPTYPE>
struct cal_force_ew_sincos_op<FPTYPE, base_device::DEVICE_GPU>
{
    void operator()(const base_device::DEVICE_GPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force);
};

#endif // __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM

} // namespace hamilt

#endif // FORCE_SINCOS_OP_H 
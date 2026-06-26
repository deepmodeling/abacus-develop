#ifndef KEDF_LKT_H
#define KEDF_LKT_H
#include <cmath>
#include <cstdio>

#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

#ifdef __CUDA
#include <cufft.h>
#endif

/**
 * @brief A class which calculates the kinetic energy, potential, and stress with Luo-Karasiev-Trickey (LKT) KEDF.
 * See Luo K, Karasiev V V, Trickey S B. Physical Review B, 2018, 98(4): 041111.
 * @author sunliang on 2023-04-28
 *
 * Optimization: pre-allocated working buffers to eliminate per-call new/delete overhead;
 * added OpenMP parallelization for real-space grid loops.
 */
class KEDF_LKT
{
  public:
    KEDF_LKT()
    {
        this->stress.create(3, 3);
    }
    ~KEDF_LKT()
    {
        free_buffers();
#ifdef __CUDA
        free_gpu_buffers();
#endif
    }

    void set_para(double dV, double lkt_a);

    /// Allocate persistent working buffers (called once in init)
    void init_buffers(const int nrxx);

    /// Free working buffers
    void free_buffers();

    double get_energy(const double* const* prho, ModulePW::PW_Basis* pw_rho);
    double get_energy_density(const double* const* prho, int is, int ir, ModulePW::PW_Basis* pw_rho);
    void tau_lkt(const double* const* prho, ModulePW::PW_Basis* pw_rho, double* rtau_lkt);
    void lkt_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpotential);
    void get_stress(const double* const* prho, ModulePW::PW_Basis* pw_rho);

    double lkt_energy = 0.; // LKT energy
    ModuleBase::matrix stress;

  private:
    void nabla(const double* pinput, ModulePW::PW_Basis* pw_rho, double** routput);
    void divergence(const double* const* pinput, ModulePW::PW_Basis* pw_rho, double* routput);
    void get_as(const double* prho, const double* const* pnabla_rho, const int nrxx, double* as);

    double dV_ = 0.; // volume element = V/nxyz
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double s_coef_
        = 1.0 / (2. * std::pow(3 * std::pow(M_PI, 2.0), 1.0 / 3.0)); // coef of s, s=s_coef * |nabla rho|/rho^{4/3}
    double lkt_a_ = 1.3;

    // Pre-allocated working buffers (eliminate per-call new/delete)
    double* as_ = nullptr;             // a*s values, size nrxx
    double* nabla_rho_[3] = {nullptr, nullptr, nullptr}; // gradient components
    double* div_input_[3] = {nullptr, nullptr, nullptr}; // input for divergence()
    double* nabla_term_ = nullptr;     // divergence output, size nrxx
    int buffer_nrxx_ = 0;              // size of allocated buffers

#ifdef __CUDA
    // ── GPU acceleration ──
    /// GPU-accelerated lkt_potential (cuFFT-based gradient + divergence).
    void lkt_potential_gpu(const double* const* prho, ModulePW::PW_Basis* pw_rho,
                           ModuleBase::matrix& rpotential);

    /// Release persistent GPU buffers and cuFFT plans.
    void free_gpu_buffers();

    // Persistent GPU buffers (lazily allocated once, reused across SCF iterations)
    double* d_rho_ = nullptr;            // real-space density       (nrxx doubles)
    double* d_as_ = nullptr;             // a*s values               (nrxx doubles)
    double* d_potential_ = nullptr;      // output potential         (nrxx doubles)
    double* d_fft_save_ = nullptr;       // cuFFT complex workspace  (nrxx×2 doubles)
    double* d_fft_work_ = nullptr;       // cuFFT complex workspace  (nrxx×2 doubles)
    double* d_nabla_rho_[3] = {nullptr, nullptr, nullptr}; // gradient components (nrxx doubles each)
    double* d_div_input_[3] = {nullptr, nullptr, nullptr}; // divergence input   (alias d_nabla_rho_)
    double* d_gcar_ = nullptr;           // interleaved G vectors    (npw×3 doubles)
    double* d_energy_partial_ = nullptr; // energy partial sums      (nblocks doubles)

    cufftHandle cufft_plan_fwd_ = 0;    // cuFFT forward plan
    cufftHandle cufft_plan_bwd_ = 0;    // cuFFT backward plan

    bool gpu_allocated_ = false;
#endif
};
#endif

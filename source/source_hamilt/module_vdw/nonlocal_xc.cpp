#include "nonlocal_xc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

// vdwxc.h includes <complex.h> inside an extern "C" block. Include the C++ header first so that
// libc++ does not parse its templates with C linkage, then drop the I macro.
#include <complex.h>
#include <complex>
#include <vdwxc.h>
#ifdef __MPI
#include <vdwxc_mpi.h>
#endif
#undef I

#include <algorithm>
#include <cassert>
#include <vector>

namespace vdw
{
namespace
{

/// Cartesian components of the stress tensor returned by vdwxc_stress, stored as [3 * a + b].
const int stress_components = 9;

/// Spin channels of chr.rho that add up to the valence density (for nspin = 4, rho[0] is the total density).
int density_channels(const Charge &chr) { return (chr.nspin == 2) ? 2 : 1; }

/// z planes [begin, end) of the libvdwxc slab of process rank: FFTW-MPI blocks of ceil(nz / nproc) planes.
void libvdwxc_slab(
    const int nz,
    const int nproc,
    const int rank,
    int &begin,
    int &end)
{
    const int block = (nz + nproc - 1) / nproc;
    const int first_plane = rank * block;
    begin = std::min(first_plane, nz);
    const int block_end = begin + block;
    end = std::min(block_end, nz);
}

/**
 * @brief Moves a field between the ABACUS distribution of pw and the libvdwxc slabs of the pool.
 *
 * libvdwxc gets the lattice vectors in the cyclic order (a3, a1, a2), so that the first axis of its grid is the
 * z axis that ABACUS distributes. It splits that axis into blocks of ceil(nz/P) planes, the FFTW-MPI default;
 * vdwxc_init_mpi does not report the local slab, so the same rule is applied here. The slab of a process is
 * stored as [z][x][y], the transpose of the ABACUS layout [x][y][z], and whole z planes are moved between the
 * two distributions by all-to-all exchanges.
 */
class SlabExchange
{
  public:
    explicit SlabExchange(const ModulePW::PW_Basis &pw) : pw_(pw)
    {
        const int nproc = pw.poolnproc;
        int slab_start = 0;
        int slab_end = 0;
        libvdwxc_slab(pw.nz, nproc, pw.poolrank, slab_start, slab_end);
        const int plane_end = pw.startz_current + pw.nplane;
        send_count_.assign(nproc, 0);
        send_offset_.assign(nproc, 0);
        recv_count_.assign(nproc, 0);
        recv_offset_.assign(nproc, 0);
        for (int ip = 0; nproc > 1 && ip < nproc; ++ip)
        {
            // ABACUS planes of this process that lie in the slab of ip
            int slab_start_of_ip = 0;
            int slab_end_of_ip = 0;
            libvdwxc_slab(pw.nz, nproc, ip, slab_start_of_ip, slab_end_of_ip);
            int low = std::max(pw.startz_current, slab_start_of_ip);
            int high = std::min(plane_end, slab_end_of_ip);
            if (high > low)
            {
                send_count_[ip] = (high - low) * pw.nxy;
                send_offset_[ip] = (low - pw.startz_current) * pw.nxy;
            }
            // ABACUS planes of ip that lie in the slab of this process
            const int plane_end_of_ip = pw.startz[ip] + pw.numz[ip];
            low = std::max(pw.startz[ip], slab_start);
            high = std::min(plane_end_of_ip, slab_end);
            if (high > low)
            {
                recv_count_[ip] = (high - low) * pw.nxy;
                recv_offset_[ip] = (low - slab_start) * pw.nxy;
            }
        }
        // libvdwxc expects valid pointers also for an empty slab
        const int slab_points = (slab_end - slab_start) * pw.nxy;
        slab_size_ = std::max(1, slab_points);
        // the transposed ABACUS planes are only needed for the exchange between processes
        const int planes_size = (nproc > 1) ? std::max(1, pw.nrxx) : 0;
        planes_.resize(planes_size);
    }

    /// Number of doubles of the libvdwxc slab of this process.
    int slab_size() const { return slab_size_; }

    /// field on the nrxx points of this process to the libvdwxc slab
    void to_slab(const std::vector<double> &field, std::vector<double> &slab)
    {
        std::vector<double> &planes = (pw_.poolnproc > 1) ? planes_ : slab;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ixy = 0; ixy < pw_.nxy; ++ixy)
        {
            for (int iz = 0; iz < pw_.nplane; ++iz)
            {
                planes[iz * pw_.nxy + ixy] = field[ixy * pw_.nplane + iz];
            }
        }
#ifdef __MPI
        if (pw_.poolnproc > 1)
        {
            Parallel_Common::alltoallv_data(
                planes_.data(), send_count_.data(), send_offset_.data(), slab.data(), recv_count_.data(),
                recv_offset_.data(), pw_.pool_world);
        }
#endif
    }

    /// libvdwxc slab to field on the nrxx points of this process
    void to_grid(const std::vector<double> &slab, std::vector<double> &field)
    {
#ifdef __MPI
        if (pw_.poolnproc > 1)
        {
            Parallel_Common::alltoallv_data(
                slab.data(), recv_count_.data(), recv_offset_.data(), planes_.data(), send_count_.data(),
                send_offset_.data(), pw_.pool_world);
        }
#endif
        const std::vector<double> &planes = (pw_.poolnproc > 1) ? planes_ : slab;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ixy = 0; ixy < pw_.nxy; ++ixy)
        {
            for (int iz = 0; iz < pw_.nplane; ++iz)
            {
                field[ixy * pw_.nplane + iz] = planes[iz * pw_.nxy + ixy];
            }
        }
    }

  private:
    const ModulePW::PW_Basis &pw_;
    // planes moved to (send) and from (recv) each process, in doubles
    std::vector<int> send_count_;
    std::vector<int> send_offset_;
    std::vector<int> recv_count_;
    std::vector<int> recv_offset_;
    int slab_size_ = 1;
    // the ABACUS planes of this process as [z][xy]
    std::vector<double> planes_;
};

/**
 * @brief rVV10 energy and derivatives of an ABACUS-distributed density, from libvdwxc.
 *
 * The libvdwxc object lives for this one evaluation, like the Libxc functionals of an XC evaluation: its
 * kernel table is static and its FFT plans are FFTW_ESTIMATE plans, so making it is cheap next to the
 * evaluation, and nothing is left to release at the end of the run.
 *
 * @param pw          grid and pool of the density (chr.rhopw)
 * @param parameters  the rVV10 parameters b and C
 * @param rho         total density on the nrxx points of this process
 * @param sigma       |grad rho|^2 on the same points
 * @param dedn        dE/dn on the same points, or nullptr when it is not needed; the choice is collective
 * @param dedsigma    dE/dsigma on the same points
 * @param stress      with stress != nullptr also the explicit strain derivative, 9 Cartesian components
 *                    (Ha, not divided by the volume)
 * @return the energy in Ha, reduced over the pool
 */
double evaluate_libvdwxc(
    const ModulePW::PW_Basis &pw,
    const XCFunctionalParameters &parameters,
    const std::vector<double> &rho,
    const std::vector<double> &sigma,
    std::vector<double> *dedn,
    std::vector<double> &dedsigma,
    double *stress)
{
    ModuleBase::timer::start("vdw", "evaluate_libvdwxc");
    vdwxc_data data = vdwxc_new(VDWXC_RVV10);
    if (vdwxc_set_rvv10_parameters(data, parameters.rvv10_b, parameters.rvv10_c) != 0)
    {
        ModuleBase::WARNING_QUIT("vdw::evaluate_libvdwxc", "libvdwxc rejected rvv10_b or rvv10_c");
    }
    // libvdwxc integrates with the signed cell volume, and inverting a left-handed cell leaves rVV10 unchanged
    ModuleBase::Matrix3 a = pw.latvec * pw.lat0;
    if (a.Det() < 0.0)
    {
        a = a * -1.0;
    }
    // the cell must be known before the FFT plans are made
    vdwxc_set_unit_cell(data, pw.nz, pw.nx, pw.ny, a.e31, a.e32, a.e33, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23);
    const bool distributed = pw.poolnproc > 1;
#ifdef __MPI
    if (distributed)
    {
        vdwxc_init_mpi(data, pw.pool_world);
    }
#endif
    if (!distributed)
    {
        vdwxc_init_serial(data);
    }

    SlabExchange exchange(pw);
    const int slab_size = exchange.slab_size();
    std::vector<double> slab_rho(slab_size);
    std::vector<double> slab_sigma(slab_size);
    std::vector<double> slab_dedn(slab_size);
    std::vector<double> slab_dedsigma(slab_size);
    exchange.to_slab(rho, slab_rho);
    exchange.to_slab(sigma, slab_sigma);
    double energy = 0.0;
    if (stress == nullptr)
    {
        energy = vdwxc_calculate(data, slab_rho.data(), slab_sigma.data(), slab_dedn.data(), slab_dedsigma.data());
    }
    else
    {
        energy = vdwxc_stress(data, stress, slab_rho.data(), slab_sigma.data(), slab_dedn.data(), slab_dedsigma.data());
    }
    vdwxc_finalize(&data);
    if (dedn != nullptr)
    {
        exchange.to_grid(slab_dedn, *dedn);
    }
    exchange.to_grid(slab_dedsigma, dedsigma);
#ifdef __MPI
    // the distributed backend returns the contribution of this process
    if (distributed)
    {
        Parallel_Common::reduce_data(&energy, 1, pw.pool_world);
        if (stress != nullptr)
        {
            Parallel_Common::reduce_data(stress, stress_components, pw.pool_world);
        }
    }
#endif
    ModuleBase::timer::end("vdw", "evaluate_libvdwxc");
    return energy;
}

/// rVV10 energy (Ha) of the total density of chr (valence plus core), with the plane-wave gradient of that
/// density (a copy of total_gradient when that is given), dE/dn (unless dedn is nullptr), dE/dsigma for
/// sigma = |grad n|^2 and, with stress != nullptr, the explicit strain derivative.
double evaluate_total_density(
    const XCFunctionalParameters &parameters,
    const Charge &chr,
    const std::vector<ModuleBase::Vector3<double>> *total_gradient,
    std::vector<ModuleBase::Vector3<double>> &grad,
    std::vector<double> *dedn,
    std::vector<double> &dedsigma,
    double *stress)
{
    const ModulePW::PW_Basis &pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    const int nspin0 = density_channels(chr);
    const double *core_begin = chr.rho_core;
    const double *core_end = core_begin + nrxx;
    std::vector<double> rho(core_begin, core_end);
    for (int is = 0; is < nspin0; ++is)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rho[ir] += chr.rho[is][ir];
        }
    }
    // the branch is collective (the transforms are): the caller decides it for the whole pool
    if (total_gradient != nullptr)
    {
        assert(total_gradient->size() == static_cast<std::size_t>(nrxx));
        grad = *total_gradient;
    }
    else
    {
        std::vector<std::complex<double>> rhog(pw.npw);
        pw.real2recip(rho.data(), rhog.data());
        grad.resize(nrxx);
        XC_Functional::grad_rho(rhog.data(), grad.data(), &pw, pw.tpiba);
    }
    std::vector<double> sigma(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }
    if (dedn != nullptr)
    {
        dedn->resize(nrxx);
    }
    dedsigma.resize(nrxx);
    return evaluate_libvdwxc(pw, parameters, rho, sigma, dedn, dedsigma, stress);
}
} // namespace

void add_rvv10(
    const XCFunctionalParameters &parameters,
    const Charge &chr,
    const std::vector<ModuleBase::Vector3<double>> *total_gradient,
    double &etxc,
    double &vtxc,
    ModuleBase::matrix &v)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10");
    ModuleBase::timer::start("vdw", "add_rvv10");
    const ModulePW::PW_Basis &pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> dedn;
    std::vector<double> dedsigma;
    const double energy = evaluate_total_density(parameters, chr, total_gradient, grad, &dedn, dedsigma, nullptr);
    // v = dE/dn - div(2 dE/dsigma grad n), Hartree to Ry by e2
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &pw, pw.tpiba);
    // for nspin = 4 the potential of the total density goes to v(0) only
    const int nspin0 = density_channels(chr);
    // this reduction and the one of the stress stay serial, so the result does not depend on the thread count
    double vtxc_nl = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double v_nl = ModuleBase::e2 * (dedn[ir] - divergence[ir]);
        for (int is = 0; is < nspin0; ++is)
        {
            v(is, ir) += v_nl;
            vtxc_nl += v_nl * chr.rho[is][ir];
        }
    }
#ifdef __MPI
    Parallel_Common::reduce_data(&vtxc_nl, 1, pw.pool_world);
#endif
    etxc += ModuleBase::e2 * energy;
    vtxc += vtxc_nl * pw.omega / pw.nxyz;
    ModuleBase::timer::end("vdw", "add_rvv10");
}

void add_rvv10_stress(const XCFunctionalParameters &parameters, const Charge &chr, ModuleBase::matrix &sigma)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10_stress");
    ModuleBase::timer::start("vdw", "add_rvv10_stress");
    const ModulePW::PW_Basis &pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> dedsigma;
    double explicit_stress[stress_components] = {0.0};
    // the stress needs dE/dsigma only
    const double energy = evaluate_total_density(parameters, chr, nullptr, grad, nullptr, dedsigma, explicit_stress);
    // integral of dE/dsigma d_a n d_b n, in units of the volume element
    double gradient_term[stress_components] = {0.0};
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double g[3] = {grad[ir].x, grad[ir].y, grad[ir].z};
        for (int a = 0; a < 3; ++a)
        {
            for (int b = 0; b < 3; ++b)
            {
                gradient_term[3 * a + b] += dedsigma[ir] * g[a] * g[b];
            }
        }
    }
#ifdef __MPI
    Parallel_Common::reduce_data(gradient_term, stress_components, pw.pool_world);
#endif
    // The full strain derivative is T - delta_ab Int(n dE/dn + 2 sigma dE/dsigma) - 2 Int dE/dsigma d_a n d_b n,
    // T the explicit part from libvdwxc. The XC diagonal -(etxc - vtxc)/omega in PW (in LCAO -etxc/omega plus the
    // Pulay term of the potential) already holds -e2 (E - Int(n dE/dn + 2 sigma dE/dsigma)) / omega; the rest is below.
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            const double diagonal = (a == b) ? energy : 0.0;
            sigma(a, b) += -ModuleBase::e2 / pw.omega * (explicit_stress[3 * a + b] - diagonal) +
                           2.0 * ModuleBase::e2 * gradient_term[3 * a + b] / pw.nxyz;
        }
    }
    ModuleBase::timer::end("vdw", "add_rvv10_stress");
}

} // namespace vdw

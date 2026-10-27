#include "nonlocal_xc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

#ifdef __MPI
// vdwxc_mpi.h includes <mpi.h> inside an extern "C" block, so include it first here.
#include <mpi.h>
#endif
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
#include <vector>

namespace vdw
{
namespace
{
/// Cartesian components of the stress tensor returned by vdwxc_stress, stored as [3 * a + b]
const int STRESS_COMPONENTS = 9;

/**
 * libvdwxc object for one grid and pool, kept between calls.
 *
 * libvdwxc distributes the first axis of its grid in blocks of ceil(n/P) planes, the FFTW-MPI default
 * (vdwxc_init_mpi does not report its local slab, so plan_exchange recomputes it).
 * It is given the lattice vectors in the cyclic order (a3, a1, a2), so that this axis is the z axis
 * that ABACUS distributes and the slab of a process is stored as [z][x][y], the transpose of the
 * ABACUS layout [x][y][z]. Whole z planes are moved between the two distributions with all-to-all exchanges.
 */
class Rvv10Workspace
{
  public:
    /// libvdwxc energy (Ha) and derivatives of the ABACUS-distributed fields rho and sigma; with
    /// stress != nullptr also the explicit strain derivative (Ha, Cartesian, not divided by the volume).
    double evaluate(const ModulePW::PW_Basis& pw,
                    const double b,
                    const double c,
                    const std::vector<double>& rho,
                    const std::vector<double>& sigma,
                    std::vector<double>& dedn,
                    std::vector<double>& dedsigma,
                    double* stress);

  private:
    void prepare(const ModulePW::PW_Basis& pw, const double b, const double c);
    void plan_exchange(const ModulePW::PW_Basis& pw);
    void to_libvdwxc(const ModulePW::PW_Basis& pw, const std::vector<double>& field, std::vector<double>& slab);
    void from_libvdwxc(const ModulePW::PW_Basis& pw, const std::vector<double>& slab, std::vector<double>& field);

    vdwxc_data data_ = nullptr;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    int nproc_ = 0;
#ifdef __MPI
    MPI_Comm comm_ = MPI_COMM_NULL;
#endif
    // planes moved to (send) and from (recv) each process, in doubles
    std::vector<int> send_count_;
    std::vector<int> send_offset_;
    std::vector<int> recv_count_;
    std::vector<int> recv_offset_;
    // the ABACUS planes of this process as [z][xy], and the libvdwxc slab
    std::vector<double> planes_;
    std::vector<double> rho_;
    std::vector<double> sigma_;
    std::vector<double> dedn_;
    std::vector<double> dedsigma_;
};

void Rvv10Workspace::prepare(const ModulePW::PW_Basis& pw, const double b, const double c)
{
    bool same_grid = data_ != nullptr && pw.nx == nx_ && pw.ny == ny_ && pw.nz == nz_ && pw.poolnproc == nproc_;
#ifdef __MPI
    same_grid = same_grid && (pw.poolnproc == 1 || pw.pool_world == comm_);
#endif
    if (!same_grid)
    {
        if (data_ != nullptr)
        {
            vdwxc_finalize(&data_);
        }
        data_ = vdwxc_new(VDWXC_RVV10);
        nx_ = pw.nx;
        ny_ = pw.ny;
        nz_ = pw.nz;
        nproc_ = pw.poolnproc;
    }
    // The cell is set on every call: a cell relaxation keeps the FFT grid but changes the cell. libvdwxc
    // integrates with the signed cell volume, and inverting a left-handed cell leaves rVV10 unchanged.
    ModuleBase::Matrix3 a = pw.latvec * pw.lat0;
    if (a.Det() < 0.0)
    {
        a = a * -1.0;
    }
    vdwxc_set_unit_cell(data_, pw.nz, pw.nx, pw.ny, a.e31, a.e32, a.e33, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23);
    if (!same_grid)
    {
#ifdef __MPI
        comm_ = pw.pool_world;
        if (pw.poolnproc > 1)
        {
            vdwxc_init_mpi(data_, comm_);
        }
        else
        {
            vdwxc_init_serial(data_);
        }
#else
        vdwxc_init_serial(data_);
#endif
    }
    // The libvdwxc slabs depend only on nz and the number of processes, but the ABACUS planes depend on
    // the basis, so the exchange is planned on every call.
    this->plan_exchange(pw);
    if (vdwxc_set_rvv10_parameters(data_, b, c) != 0)
    {
        ModuleBase::WARNING_QUIT("vdw::add_rvv10", "libvdwxc rejected vdw_rvv10_b or vdw_rvv10_c");
    }
}

void Rvv10Workspace::plan_exchange(const ModulePW::PW_Basis& pw)
{
    const int nproc = pw.poolnproc;
    const int block = (pw.nz + nproc - 1) / nproc;
    const int first_slab_plane = pw.poolrank * block;
    const int slab_start = std::min(first_slab_plane, pw.nz);
    const int planes_after_start = pw.nz - slab_start;
    const int nslab = std::min(block, planes_after_start);
    const int slab_end = slab_start + nslab;
    const int plane_end = pw.startz_current + pw.nplane;
    send_count_.assign(nproc, 0);
    send_offset_.assign(nproc, 0);
    recv_count_.assign(nproc, 0);
    recv_offset_.assign(nproc, 0);
    for (int ip = 0; nproc > 1 && ip < nproc; ++ip)
    {
        // ABACUS planes of this process that lie in the slab of ip
        const int first_plane_of_ip = ip * block;
        const int begin = std::min(first_plane_of_ip, pw.nz);
        const int block_end = begin + block;
        const int end = std::min(block_end, pw.nz);
        int low = std::max(pw.startz_current, begin);
        int high = std::min(plane_end, end);
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
    const int slab_points = nslab * pw.nxy;
    const int slab_size = std::max(1, slab_points);
    // the transposed ABACUS planes are only needed for the exchange between processes
    const int planes_size = (nproc > 1) ? std::max(1, pw.nrxx) : 0;
    planes_.resize(planes_size);
    rho_.resize(slab_size);
    sigma_.resize(slab_size);
    dedn_.resize(slab_size);
    dedsigma_.resize(slab_size);
}

void Rvv10Workspace::to_libvdwxc(const ModulePW::PW_Basis& pw,
                                 const std::vector<double>& field,
                                 std::vector<double>& slab)
{
    std::vector<double>& planes = (pw.poolnproc > 1) ? planes_ : slab;
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            planes[iz * pw.nxy + ixy] = field[ixy * pw.nplane + iz];
        }
    }
#ifdef __MPI
    if (pw.poolnproc > 1)
    {
        Parallel_Common::alltoallv_data(planes_.data(),
                                        send_count_.data(),
                                        send_offset_.data(),
                                        slab.data(),
                                        recv_count_.data(),
                                        recv_offset_.data(),
                                        pw.pool_world);
    }
#endif
}

void Rvv10Workspace::from_libvdwxc(const ModulePW::PW_Basis& pw,
                                   const std::vector<double>& slab,
                                   std::vector<double>& field)
{
#ifdef __MPI
    if (pw.poolnproc > 1)
    {
        Parallel_Common::alltoallv_data(slab.data(),
                                        recv_count_.data(),
                                        recv_offset_.data(),
                                        planes_.data(),
                                        send_count_.data(),
                                        send_offset_.data(),
                                        pw.pool_world);
    }
#endif
    const std::vector<double>& planes = (pw.poolnproc > 1) ? planes_ : slab;
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            field[ixy * pw.nplane + iz] = planes[iz * pw.nxy + ixy];
        }
    }
}

double Rvv10Workspace::evaluate(const ModulePW::PW_Basis& pw,
                                const double b,
                                const double c,
                                const std::vector<double>& rho,
                                const std::vector<double>& sigma,
                                std::vector<double>& dedn,
                                std::vector<double>& dedsigma,
                                double* stress)
{
    this->prepare(pw, b, c);
    this->to_libvdwxc(pw, rho, rho_);
    this->to_libvdwxc(pw, sigma, sigma_);
    double energy = 0.0;
    if (stress == nullptr)
    {
        energy = vdwxc_calculate(data_, rho_.data(), sigma_.data(), dedn_.data(), dedsigma_.data());
    }
    else
    {
        energy = vdwxc_stress(data_, stress, rho_.data(), sigma_.data(), dedn_.data(), dedsigma_.data());
    }
    this->from_libvdwxc(pw, dedn_, dedn);
    this->from_libvdwxc(pw, dedsigma_, dedsigma);
#ifdef __MPI
    // the distributed backend returns the contribution of this process
    if (pw.poolnproc > 1)
    {
        Parallel_Common::reduce_data(&energy, 1, pw.pool_world);
        if (stress != nullptr)
        {
            Parallel_Common::reduce_data(stress, STRESS_COMPONENTS, pw.pool_world);
        }
    }
#endif
    return energy;
}

/// The workspace is never destroyed: a static object would be destroyed after MPI_Finalize and
/// fftw_cleanup_threads, which libvdwxc's FFT plans must not outlive.
Rvv10Workspace& workspace()
{
    static Rvv10Workspace* const instance = new Rvv10Workspace;
    return *instance;
}

/// Total density of chr (valence plus core), its plane-wave gradient and sigma = |grad n|^2.
void total_density(const Charge& chr,
                   std::vector<double>& rho,
                   std::vector<ModuleBase::Vector3<double>>& grad,
                   std::vector<double>& sigma)
{
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    const double* core_begin = chr.rho_core;
    const double* core_end = core_begin + nrxx;
    // for nspin = 4, rho[0] is the total valence density
    const int nspin0 = (chr.nspin == 2) ? 2 : 1;
    rho.assign(core_begin, core_end);
    for (int is = 0; is < nspin0; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rho[ir] += chr.rho[is][ir];
        }
    }
    std::vector<std::complex<double>> rhog(pw.npw);
    pw.real2recip(rho.data(), rhog.data());
    grad.resize(nrxx);
    XC_Functional::grad_rho(rhog.data(), grad.data(), &pw, pw.tpiba);
    sigma.resize(nrxx);
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }
}
} // namespace

void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10");
    ModuleBase::timer::start("vdw", "add_rvv10");
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<double> rho;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> sigma;
    total_density(chr, rho, grad, sigma);
    std::vector<double> dedn(nrxx);
    std::vector<double> dedsigma(nrxx);
    const double energy
        = workspace().evaluate(pw, parameters.rvv10_b, parameters.rvv10_c, rho, sigma, dedn, dedsigma, nullptr);

    // v = dE/dn - div(2 dE/dsigma grad n), converted from Hartree to Ry by e2
    for (int ir = 0; ir < nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &pw, pw.tpiba);
    // for nspin = 4 the potential of the total density goes to v(0) only
    const int nspin0 = (chr.nspin == 2) ? 2 : 1;
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

void add_rvv10_stress(const XCFunctionalParameters& parameters, const Charge& chr, ModuleBase::matrix& sigma)
{
    if (!parameters.rvv10)
    {
        return;
    }
    ModuleBase::TITLE("vdw", "add_rvv10_stress");
    ModuleBase::timer::start("vdw", "add_rvv10_stress");
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    std::vector<double> rho;
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> density_sigma;
    total_density(chr, rho, grad, density_sigma);
    std::vector<double> dedn(nrxx);
    std::vector<double> dedsigma(nrxx);
    double explicit_stress[STRESS_COMPONENTS] = {0.0};
    const double energy = workspace().evaluate(pw,
                                               parameters.rvv10_b,
                                               parameters.rvv10_c,
                                               rho,
                                               density_sigma,
                                               dedn,
                                               dedsigma,
                                               explicit_stress);

    // integral of dE/dsigma d_a n d_b n, in units of the volume element
    double gradient_term[STRESS_COMPONENTS] = {0.0};
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
    Parallel_Common::reduce_data(gradient_term, STRESS_COMPONENTS, pw.pool_world);
#endif

    // The full strain derivative is T - delta_ab Int(n dE/dn + 2 sigma dE/dsigma) - 2 Int dE/dsigma d_a n d_b n,
    // with T the explicit part from libvdwxc. The XC diagonal -(etxc - vtxc)/omega in PW (in LCAO, -etxc/omega
    // and the Pulay term of the effective potential) already holds -e2 (E - Int(n dE/dn + 2 sigma dE/dsigma)) / omega,
    // which leaves the terms below (ABACUS sign, Ry).
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            const double diagonal = (a == b) ? energy : 0.0;
            sigma(a, b) += -ModuleBase::e2 / pw.omega * (explicit_stress[3 * a + b] - diagonal)
                           + 2.0 * ModuleBase::e2 * gradient_term[3 * a + b] / pw.nxyz;
        }
    }
    ModuleBase::timer::end("vdw", "add_rvv10_stress");
}

} // namespace vdw

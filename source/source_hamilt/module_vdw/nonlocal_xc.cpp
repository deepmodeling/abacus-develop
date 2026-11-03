#include "nonlocal_xc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/memory_recorder.h"
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
#include <cassert>
#include <vector>

namespace vdw
{
namespace
{
/// Cartesian components of the stress tensor returned by vdwxc_stress, stored as [3 * a + b]
const int stress_components = 9;

/// Spin channels of chr.rho that add up to the valence density (for nspin = 4, rho[0] is the total density).
int density_channels(const Charge& chr)
{
    return (chr.nspin == 2) ? 2 : 1;
}

/// z planes [begin, end) of the libvdwxc slab of process rank: FFTW-MPI blocks of ceil(nz / nproc) planes.
void libvdwxc_slab(const int nz, const int nproc, const int rank, int& begin, int& end)
{
    const int block = (nz + nproc - 1) / nproc;
    const int first_plane = rank * block;
    begin = std::min(first_plane, nz);
    const int block_end = begin + block;
    end = std::min(block_end, nz);
}

/**
 * libvdwxc object for one FFT grid and pool, kept between calls and rebuilt when either changes.
 * libvdwxc gets the lattice vectors in the cyclic order (a3, a1, a2), so that the first axis of its grid
 * is the z axis that ABACUS distributes. It splits that axis into blocks of ceil(nz/P) planes, the FFTW-MPI
 * default; vdwxc_init_mpi does not report the local slab, so libvdwxc_slab computes it the same way. The
 * slab of a process is stored as [z][x][y], the transpose of the ABACUS layout [x][y][z], and whole z
 * planes are moved between the two distributions by all-to-all exchanges.
 * No user-defined destructor on purpose: the buffers are freed at exit, but data_ is never finalized,
 * because the single instance (see workspace()) is destroyed after MPI_Finalize and fftw_cleanup_threads,
 * after which the FFTW-MPI plans of data_ can no longer be destroyed safely.
 */
class Rvv10Workspace
{
  public:
    Rvv10Workspace() = default;
    Rvv10Workspace(const Rvv10Workspace&) = delete;
    Rvv10Workspace& operator=(const Rvv10Workspace&) = delete;

    /// libvdwxc energy (Ha) and derivatives of the ABACUS-distributed fields rho and sigma. dedn may be
    /// nullptr when dE/dn is not needed; with stress != nullptr also the explicit strain derivative
    /// (Ha, Cartesian, not divided by the volume).
    double evaluate(const ModulePW::PW_Basis& pw,
                    const XCFunctionalParameters& parameters,
                    const std::vector<double>& rho,
                    const std::vector<double>& sigma,
                    std::vector<double>* dedn,
                    std::vector<double>& dedsigma,
                    double* stress);

  private:
    bool same_grid(const ModulePW::PW_Basis& pw) const;
    void set_cell(const ModulePW::PW_Basis& pw);
    void rebuild(const ModulePW::PW_Basis& pw);
    void prepare(const ModulePW::PW_Basis& pw, const XCFunctionalParameters& parameters);
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

bool Rvv10Workspace::same_grid(const ModulePW::PW_Basis& pw) const
{
    bool same = data_ != nullptr && pw.nx == nx_ && pw.ny == ny_ && pw.nz == nz_ && pw.poolnproc == nproc_;
#ifdef __MPI
    same = same && (pw.poolnproc == 1 || pw.pool_world == comm_);
#endif
    return same;
}

/// The cell is set on every call: a cell relaxation keeps the FFT grid but changes the cell. libvdwxc
/// integrates with the signed cell volume, and inverting a left-handed cell leaves rVV10 unchanged.
void Rvv10Workspace::set_cell(const ModulePW::PW_Basis& pw)
{
    ModuleBase::Matrix3 a = pw.latvec * pw.lat0;
    if (a.Det() < 0.0)
    {
        a = a * -1.0;
    }
    vdwxc_set_unit_cell(data_, pw.nz, pw.nx, pw.ny, a.e31, a.e32, a.e33, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23);
}

void Rvv10Workspace::rebuild(const ModulePW::PW_Basis& pw)
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
    // the cell must be known before the FFT plans are made
    this->set_cell(pw);
#ifdef __MPI
    comm_ = pw.pool_world;
    if (pw.poolnproc > 1)
    {
        vdwxc_init_mpi(data_, comm_);
        return;
    }
#endif
    vdwxc_init_serial(data_);
}

void Rvv10Workspace::prepare(const ModulePW::PW_Basis& pw, const XCFunctionalParameters& parameters)
{
    if (this->same_grid(pw))
    {
        this->set_cell(pw);
    }
    else
    {
        this->rebuild(pw);
    }
    // the libvdwxc slabs depend only on nz and the processes, but the ABACUS planes also on the basis, so
    // the exchange is planned on every call (counts and resizes only)
    this->plan_exchange(pw);
    if (vdwxc_set_rvv10_parameters(data_, parameters.rvv10_b, parameters.rvv10_c) != 0)
    {
        ModuleBase::WARNING_QUIT("Rvv10Workspace::prepare", "libvdwxc rejected vdw_rvv10_b or vdw_rvv10_c");
    }
}

void Rvv10Workspace::plan_exchange(const ModulePW::PW_Basis& pw)
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
        int begin = 0;
        int end = 0;
        libvdwxc_slab(pw.nz, nproc, ip, begin, end);
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
    const int slab_size = std::max(1, (slab_end - slab_start) * pw.nxy);
    // the transposed ABACUS planes are only needed for the exchange between processes
    const int planes_size = (nproc > 1) ? std::max(1, pw.nrxx) : 0;
    planes_.resize(planes_size);
    rho_.resize(slab_size);
    sigma_.resize(slab_size);
    dedn_.resize(slab_size);
    dedsigma_.resize(slab_size);
    const long buffer_bytes = static_cast<long>(sizeof(double)) * (planes_size + 4L * slab_size);
    ModuleBase::Memory::record("Rvv10Workspace::buffers", buffer_bytes);
}

void Rvv10Workspace::to_libvdwxc(const ModulePW::PW_Basis& pw,
                                 const std::vector<double>& field,
                                 std::vector<double>& slab)
{
    std::vector<double>& planes = (pw.poolnproc > 1) ? planes_ : slab;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            field[ixy * pw.nplane + iz] = planes[iz * pw.nxy + ixy];
        }
    }
}

double Rvv10Workspace::evaluate(const ModulePW::PW_Basis& pw,
                                const XCFunctionalParameters& parameters,
                                const std::vector<double>& rho,
                                const std::vector<double>& sigma,
                                std::vector<double>* dedn,
                                std::vector<double>& dedsigma,
                                double* stress)
{
    ModuleBase::timer::start("Rvv10Workspace", "evaluate");
    this->prepare(pw, parameters);
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
    // the exchange is collective, and every process of the pool asks for dE/dn or none does
    if (dedn != nullptr)
    {
        this->from_libvdwxc(pw, dedn_, *dedn);
    }
    this->from_libvdwxc(pw, dedsigma_, dedsigma);
#ifdef __MPI
    if (pw.poolnproc > 1) // the distributed backend returns the contribution of this process
    {
        Parallel_Common::reduce_data(&energy, 1, pw.pool_world);
        if (stress != nullptr)
        {
            Parallel_Common::reduce_data(stress, stress_components, pw.pool_world);
        }
    }
#endif
    ModuleBase::timer::end("Rvv10Workspace", "evaluate");
    return energy;
}

Rvv10Workspace& workspace()
{
    static Rvv10Workspace instance;
    return instance;
}

/// Total density of chr (valence plus core), its plane-wave gradient (a copy of total_gradient when that
/// is given) and sigma = |grad n|^2.
void total_density(const Charge& chr,
                   const std::vector<ModuleBase::Vector3<double>>* total_gradient,
                   std::vector<double>& rho,
                   std::vector<ModuleBase::Vector3<double>>& grad,
                   std::vector<double>& sigma)
{
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    const int nspin0 = density_channels(chr);
    rho.assign(chr.rho_core, chr.rho_core + nrxx);
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
    // the branch is collective (the transforms below are): the caller decides it for the whole pool
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
    sigma.resize(nrxx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1024)
#endif
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }
}

/// rVV10 energy (Ha) of the total density of chr, with the density gradient (see total_density), dE/dn
/// (unless dedn is nullptr), dE/dsigma and, with stress != nullptr, the explicit strain derivative.
double evaluate_total_density(const XCFunctionalParameters& parameters,
                              const Charge& chr,
                              const std::vector<ModuleBase::Vector3<double>>* total_gradient,
                              std::vector<ModuleBase::Vector3<double>>& grad,
                              std::vector<double>* dedn,
                              std::vector<double>& dedsigma,
                              double* stress)
{
    std::vector<double> rho;
    std::vector<double> sigma;
    total_density(chr, total_gradient, rho, grad, sigma);
    if (dedn != nullptr)
    {
        dedn->resize(chr.nrxx);
    }
    dedsigma.resize(chr.nrxx);
    return workspace().evaluate(*chr.rhopw, parameters, rho, sigma, dedn, dedsigma, stress);
}
} // namespace

void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               double& etxc,
               double& vtxc,
               ModuleBase::matrix& v)
{
    add_rvv10(parameters, chr, nullptr, etxc, vtxc, v);
}
void add_rvv10(const XCFunctionalParameters& parameters,
               const Charge& chr,
               const std::vector<ModuleBase::Vector3<double>>* total_gradient,
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
    double vtxc_nl = 0.0; // this reduction and the one of the stress stay serial: no thread-count dependence
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
    std::vector<ModuleBase::Vector3<double>> grad;
    std::vector<double> dedsigma;
    double explicit_stress[stress_components] = {0.0};
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
            sigma(a, b) += -ModuleBase::e2 / pw.omega * (explicit_stress[3 * a + b] - diagonal)
                           + 2.0 * ModuleBase::e2 * gradient_term[3 * a + b] / pw.nxyz;
        }
    }
    ModuleBase::timer::end("vdw", "add_rvv10_stress");
}

} // namespace vdw

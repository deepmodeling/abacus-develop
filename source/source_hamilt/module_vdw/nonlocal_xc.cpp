#include "nonlocal_xc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

// vdwxc.h includes <complex.h> inside an extern "C" block. Include the C++ header first so that
// libc++ does not parse its templates with C linkage, then drop the I macro.
#include <complex.h>
#include <complex>
#include <vdwxc.h>
#undef I

#include <vector>

namespace vdw
{
namespace
{
/// Copy the z planes [startz, startz + nplane) of a block stored as [ixy][iz] into the complete grid [ixy][z].
void put_planes(const ModulePW::PW_Basis& pw,
                const int nplane,
                const int startz,
                const double* block,
                std::vector<double>& full)
{
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < nplane; ++iz)
        {
            full[ixy * pw.nz + startz + iz] = block[ixy * nplane + iz];
        }
    }
}

/// Copy the z planes [startz, startz + nplane) of the complete grid into a block stored as [ixy][iz].
void get_planes(const ModulePW::PW_Basis& pw,
                const int nplane,
                const int startz,
                const std::vector<double>& full,
                double* block)
{
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < nplane; ++iz)
        {
            block[ixy * nplane + iz] = full[ixy * pw.nz + startz + iz];
        }
    }
}

/// The complete grid of pw in libvdwxc's [x][y][z] order on the pool root; empty on the other processes.
std::vector<double> gather_grid(const ModulePW::PW_Basis& pw, const double* local)
{
    std::vector<double> full;
#ifdef __MPI
    MPI_Comm comm = pw.pool_world;
    if (pw.poolrank != 0)
    {
        Parallel_Common::send_data(local, pw.nrxx, 0, 0, comm);
        return full;
    }
#endif
    full.resize(pw.nxyz);
    put_planes(pw, pw.nplane, pw.startz_current, local, full);
#ifdef __MPI
    std::vector<double> block;
    for (int ip = 1; ip < pw.poolnproc; ++ip)
    {
        const int count = pw.numz[ip] * pw.nxy;
        block.resize(count);
        Parallel_Common::recv_data(block.data(), count, ip, 0, comm, MPI_STATUS_IGNORE);
        put_planes(pw, pw.numz[ip], pw.startz[ip], block.data(), full);
    }
#endif
    return full;
}

/// The local z planes of pw from the complete grid held by the pool root.
void scatter_grid(const ModulePW::PW_Basis& pw, const std::vector<double>& full, double* local)
{
#ifdef __MPI
    MPI_Comm comm = pw.pool_world;
    if (pw.poolrank != 0)
    {
        Parallel_Common::recv_data(local, pw.nrxx, 0, 0, comm, MPI_STATUS_IGNORE);
        return;
    }
    std::vector<double> block;
    for (int ip = 1; ip < pw.poolnproc; ++ip)
    {
        const int count = pw.numz[ip] * pw.nxy;
        block.resize(count);
        get_planes(pw, pw.numz[ip], pw.startz[ip], full, block.data());
        Parallel_Common::send_data(block.data(), count, ip, 0, comm);
    }
#endif
    get_planes(pw, pw.nplane, pw.startz_current, full, local);
}
} // namespace

void add_rvv10(const double b,
               const double c,
               const Charge& chr,
               std::tuple<double, double, ModuleBase::matrix>& etxc_vtxc_v)
{
    ModuleBase::TITLE("vdw", "add_rvv10");
    ModuleBase::timer::start("vdw", "add_rvv10");
    const ModulePW::PW_Basis& pw = *chr.rhopw;
    const int nrxx = chr.nrxx;
    // rVV10 depends on the total density only; for nspin=4, rho[0] is the total valence density.
    const int nchannel = (chr.nspin == 2) ? 2 : 1;
    const double* core_begin = chr.rho_core;
    const double* core_end = core_begin + nrxx;
    std::vector<double> rho(core_begin, core_end);
    for (int is = 0; is < nchannel; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            rho[ir] += chr.rho[is][ir];
        }
    }

    // sigma = |grad n|^2 from the plane-wave gradient
    std::vector<std::complex<double>> rhog(pw.npw);
    pw.real2recip(rho.data(), rhog.data());
    std::vector<ModuleBase::Vector3<double>> grad(nrxx);
    XC_Functional::grad_rho(rhog.data(), grad.data(), &pw, pw.tpiba);
    std::vector<double> sigma(nrxx);
    for (int ir = 0; ir < nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }

    // libvdwxc runs serially on the pool root; energy and derivatives are in Hartree units
    std::vector<double> rho_full = gather_grid(pw, rho.data());
    std::vector<double> sigma_full = gather_grid(pw, sigma.data());
    std::vector<double> dedn_full(rho_full.size());
    std::vector<double> dedsigma_full(rho_full.size());
    double energy = 0.0;
    if (pw.poolrank == 0)
    {
        // ABACUS stores the lattice vectors as the rows of latvec; libvdwxc takes them as three
        // consecutive groups of Cartesian components in Bohr and integrates with the signed cell
        // volume, so a left-handed cell is inverted, which leaves rVV10 unchanged.
        ModuleBase::Matrix3 a = pw.latvec * pw.lat0;
        if (a.Det() < 0.0)
        {
            a = a * -1.0;
        }
        vdwxc_data data = vdwxc_new(VDWXC_RVV10);
        vdwxc_set_rvv10_parameters(data, b, c);
        vdwxc_set_unit_cell(data, pw.nx, pw.ny, pw.nz, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23, a.e31, a.e32, a.e33);
        vdwxc_init_serial(data);
        energy = vdwxc_calculate(data, rho_full.data(), sigma_full.data(), dedn_full.data(), dedsigma_full.data());
        vdwxc_finalize(&data);
    }
#ifdef __MPI
    Parallel_Common::bcast_data(&energy, 1, pw.pool_world);
#endif
    std::vector<double> dedn(nrxx);
    std::vector<double> dedsigma(nrxx);
    scatter_grid(pw, dedn_full, dedn.data());
    scatter_grid(pw, dedsigma_full, dedsigma.data());

    // v = dE/dn - div(2 dE/dsigma grad n), converted from Hartree to Ry by e2
    for (int ir = 0; ir < nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &pw, pw.tpiba);
    ModuleBase::matrix& v = std::get<2>(etxc_vtxc_v);
    double vtxc = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double v_nl = ModuleBase::e2 * (dedn[ir] - divergence[ir]);
        for (int is = 0; is < nchannel; ++is)
        {
            v(is, ir) += v_nl;
            vtxc += v_nl * chr.rho[is][ir];
        }
    }
#ifdef __MPI
    Parallel_Common::reduce_data(&vtxc, 1, pw.pool_world);
#endif
    std::get<0>(etxc_vtxc_v) += ModuleBase::e2 * energy;
    std::get<1>(etxc_vtxc_v) += vtxc * pw.omega / pw.nxyz;
    ModuleBase::timer::end("vdw", "add_rvv10");
}

} // namespace vdw

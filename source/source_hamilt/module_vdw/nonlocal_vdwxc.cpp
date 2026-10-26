#include "nonlocal_vdwxc.h"

#ifdef __LIBVDWXC

// libvdwxc's public C header includes <complex.h> inside an extern "C"
// block. Pre-include the standard header so libc++ does not parse its C++
// templates with C linkage, then drop the I macro that complex.h defines.
#include <complex.h>
#include <vdwxc.h>
#ifdef I
#undef I
#endif

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace vdw
{
namespace
{
// The cell handed to libvdwxc at construction and the cell seen at each
// evaluation come from the same PW_Basis, so only rounding may differ.
const double CELL_RELATIVE_TOLERANCE = 1.0e-12;

bool serial_cpu_grid(const ModulePW::PW_Basis& pw)
{
    return !pw.gamma_only && pw.get_device() == "cpu" && pw.get_precision() == "double" && pw.poolnproc == 1
           && pw.nrxx == pw.nxyz && pw.nx > 0 && pw.ny > 0 && pw.nz > 0 && pw.nxyz > 0 && pw.npw > 0;
}

bool right_handed_cell(const ModulePW::PW_Basis& pw)
{
    const double determinant = pw.latvec.Det();
    return std::isfinite(pw.omega) && pw.omega > 0.0 && std::isfinite(determinant) && determinant > 0.0
           && std::isfinite(pw.tpiba) && pw.tpiba > 0.0;
}

void check_serial_basis(const ModulePW::PW_Basis& pw)
{
    if (!serial_cpu_grid(pw) || !right_handed_cell(pw))
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC",
                                 "libvdwxc rVV10 requires a one-process CPU double-precision full-complex PW density "
                                 "grid and a right-handed cell");
    }
}

void check_finite(const std::vector<double>& values, const char* what)
{
    for (double value: values)
    {
        if (!std::isfinite(value))
        {
            const std::string message = std::string("libvdwxc rVV10 received a nonfinite ") + what;
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", message);
        }
    }
}

std::array<double, 9> cell_vectors(const ModulePW::PW_Basis& pw)
{
    const ModuleBase::Matrix3& a = pw.latvec;
    // ABACUS stores a1, a2, a3 as rows. The C API groups C00,C10,C20 as
    // one vector, so pass the ABACUS rows consecutively in those groups.
    return {a.e11 * pw.lat0,
            a.e12 * pw.lat0,
            a.e13 * pw.lat0,
            a.e21 * pw.lat0,
            a.e22 * pw.lat0,
            a.e23 * pw.lat0,
            a.e31 * pw.lat0,
            a.e32 * pw.lat0,
            a.e33 * pw.lat0};
}

bool same_within_rounding(const double lhs, const double rhs)
{
    double scale = 1.0;
    if (std::abs(lhs) > scale)
    {
        scale = std::abs(lhs);
    }
    if (std::abs(rhs) > scale)
    {
        scale = std::abs(rhs);
    }
    return std::abs(lhs - rhs) <= CELL_RELATIVE_TOLERANCE * scale;
}

bool same_cell(const std::array<double, 9>& lhs, const std::array<double, 9>& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (!same_within_rounding(lhs[i], rhs[i]))
        {
            return false;
        }
    }
    return true;
}

/// Host FFT gradient of the density and sigma = |grad n|^2 on the real-space grid.
void density_gradient(const ModulePW::PW_Basis& pw,
                      const std::vector<double>& density,
                      std::vector<ModuleBase::Vector3<double>>& gradient,
                      std::vector<double>& sigma)
{
    std::vector<std::complex<double>> density_g(pw.npw);
    pw.real2recip(density.data(), density_g.data());
    XC_Functional::grad_rho(density_g.data(), gradient.data(), &pw, pw.tpiba);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        sigma[r] = gradient[r].x * gradient[r].x + gradient[r].y * gradient[r].y + gradient[r].z * gradient[r].z;
    }
}

/// Energy, potential and vtxc in Ry from the library's Hartree energy and derivatives.
NonlocalXCResult assemble_result(const ModulePW::PW_Basis& pw,
                                 const double energy_ha,
                                 const std::vector<double>& dedn,
                                 const std::vector<double>& dedsigma,
                                 const std::vector<ModuleBase::Vector3<double>>& gradient,
                                 const std::vector<double>& valence_density)
{
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        if (!std::isfinite(dedn[r]) || !std::isfinite(dedsigma[r]))
        {
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc returned a nonfinite rVV10 derivative");
        }
        // Chain rule for sigma = |grad n|^2: dE/d(grad n) = 2 * dE/dsigma * grad n.
        const double coefficient = 2.0 * dedsigma[r];
        field[r].x = coefficient * gradient[r].x;
        field[r].y = coefficient * gradient[r].y;
        field[r].z = coefficient * gradient[r].z;
    }
    std::vector<double> divergence(pw.nrxx, 0.0);
    XC_Functional::grad_dot(field.data(), divergence.data(), &pw, pw.tpiba);

    // libvdwxc returns Hartree quantities; e2 converts them to Ry as in PotXC.
    // v = dE/dn - div(dE/d(grad n)); vtxc contracts v with the valence density
    // only, matching the PotXC bookkeeping. The library includes the beta term.
    NonlocalXCResult result;
    result.energy = ModuleBase::e2 * energy_ha;
    result.vtxc = 0.0;
    result.potential.assign(pw.nrxx, 0.0);
    const double dv = pw.omega / pw.nxyz;
    for (int r = 0; r < pw.nrxx; ++r)
    {
        result.potential[r] = ModuleBase::e2 * (dedn[r] - divergence[r]);
        if (!std::isfinite(result.potential[r]))
        {
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc produced a nonfinite rVV10 potential");
        }
        result.vtxc += dv * valence_density[r] * result.potential[r];
    }
    if (!std::isfinite(result.vtxc))
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc produced a nonfinite rVV10 vtxc");
    }
    return result;
}
} // namespace

LibvdwxcNonlocalXC::LibvdwxcNonlocalXC(const ModulePW::PW_Basis& pw, double b, double c)
{
    check_serial_basis(pw);
    data_ = vdwxc_new(VDWXC_RVV10);
    if (data_ == nullptr)
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc failed to create an rVV10 evaluator");
    }
    if (vdwxc_set_rvv10_parameters(data_, b, c) != 0)
    {
        vdwxc_finalize(&data_);
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc rejected the rVV10 b/C parameters");
    }

    cell_ = cell_vectors(pw);
    // libvdwxc names the arguments by matrix columns (C00,C10,C20,...).
    vdwxc_set_unit_cell(data_,
                        pw.nx,
                        pw.ny,
                        pw.nz,
                        cell_[0],
                        cell_[1],
                        cell_[2],
                        cell_[3],
                        cell_[4],
                        cell_[5],
                        cell_[6],
                        cell_[7],
                        cell_[8]);
    vdwxc_init_serial(data_);
    nx_ = pw.nx;
    ny_ = pw.ny;
    nz_ = pw.nz;
    omega_ = pw.omega;
}

LibvdwxcNonlocalXC::~LibvdwxcNonlocalXC()
{
    if (data_ != nullptr)
    {
        vdwxc_finalize(&data_);
    }
}

void LibvdwxcNonlocalXC::check_matches_initialized(const ModulePW::PW_Basis& pw) const
{
    const bool same_grid = pw.nx == nx_ && pw.ny == ny_ && pw.nz == nz_;
    const bool same_volume = same_within_rounding(pw.omega, omega_);
    if (!same_grid || !same_volume || !same_cell(cell_, cell_vectors(pw)))
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC",
                                 "libvdwxc rVV10 received a density grid different from its initialized cell");
    }
}

NonlocalXCResult LibvdwxcNonlocalXC::evaluate(const ModulePW::PW_Basis& pw,
                                              const std::vector<double>& total_density,
                                              const std::vector<double>& valence_density) const
{
    check_serial_basis(pw);
    check_matches_initialized(pw);
    const std::size_t grid_size = static_cast<std::size_t>(pw.nrxx);
    if (total_density.size() != grid_size || valence_density.size() != grid_size)
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc rVV10 received densities of the wrong size");
    }
    check_finite(total_density, "total density");
    check_finite(valence_density, "valence density");

    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<double> sigma(pw.nrxx, 0.0);
    density_gradient(pw, total_density, gradient, sigma);

    // vdwxc_calculate takes a non-const density pointer; keep the caller's array const.
    std::vector<double> density = total_density;
    std::vector<double> dedn(pw.nrxx, 0.0);
    std::vector<double> dedsigma(pw.nrxx, 0.0);
    const double energy_ha = vdwxc_calculate(data_, density.data(), sigma.data(), dedn.data(), dedsigma.data());
    if (!std::isfinite(energy_ha))
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc returned a nonfinite rVV10 energy");
    }
    return assemble_result(pw, energy_ha, dedn, dedsigma, gradient, valence_density);
}

} // namespace vdw

#endif // __LIBVDWXC

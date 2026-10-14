#include "nonlocal_xc_libvdwxc.h"

#ifdef __LIBVDWXC

#include "source_basis/module_pw/pw_basis.h"
#include "source_base/tool_quit.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <vector>

namespace vdw
{
namespace
{
void check_serial_basis(const ModulePW::PW_Basis& pw)
{
    const double determinant = pw.latvec.Det();
    if (pw.gamma_only || pw.get_device() != "cpu" || pw.get_precision() != "double" || pw.poolnproc != 1
        || pw.nrxx != pw.nxyz || pw.nx <= 0 || pw.ny <= 0 || pw.nz <= 0 || pw.nxyz <= 0 || pw.npw <= 0
        || !std::isfinite(pw.omega) || pw.omega <= 0.0 || !std::isfinite(determinant) || determinant <= 0.0
        || !std::isfinite(pw.tpiba) || pw.tpiba <= 0.0)
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC",
            "libvdwxc rVV10 requires one-process CPU/double full-complex PW density grid and a right-handed cell");
    }
}

std::array<double, 9> cell_vectors(const ModulePW::PW_Basis& pw)
{
    const auto& a = pw.latvec;
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

bool same_cell(const std::array<double, 9>& lhs, const std::array<double, 9>& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        const double lhs_abs = std::abs(lhs[i]);
        const double rhs_abs = std::abs(rhs[i]);
        double scale = 1.0;
        if (lhs_abs > scale)
        {
            scale = lhs_abs;
        }
        if (rhs_abs > scale)
        {
            scale = rhs_abs;
        }
        if (std::abs(lhs[i] - rhs[i]) > 1.0e-12 * scale)
        {
            return false;
        }
    }
    return true;
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

NonlocalXCResult LibvdwxcNonlocalXC::evaluate(const ModulePW::PW_Basis& pw,
                                              const std::vector<double>& total_density,
                                              const std::vector<double>& valence_density) const
{
    check_serial_basis(pw);
    const std::array<double, 9> current_cell = cell_vectors(pw);
    const double omega_abs = std::abs(pw.omega);
    const double initialized_omega_abs = std::abs(omega_);
    double volume_scale = 1.0;
    if (omega_abs > volume_scale)
    {
        volume_scale = omega_abs;
    }
    if (initialized_omega_abs > volume_scale)
    {
        volume_scale = initialized_omega_abs;
    }
    const bool same_grid = pw.nx == nx_ && pw.ny == ny_ && pw.nz == nz_;
    const bool same_volume = std::abs(pw.omega - omega_) <= 1.0e-12 * volume_scale;
    if (!same_grid || !same_cell(cell_, current_cell) || !same_volume
        || total_density.size() != static_cast<std::size_t>(pw.nrxx)
        || valence_density.size() != total_density.size())
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC",
                                "libvdwxc rVV10 received a density grid different from its initialized cell");
    }

    using Complex = std::complex<double>;
    std::vector<Complex> density_g(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<double> sigma(pw.nrxx, 0.0);
    std::vector<double> dedn(pw.nrxx, 0.0);
    std::vector<double> dedsigma(pw.nrxx, 0.0);
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    std::vector<double> divergence(pw.nrxx, 0.0);

    for (double value: total_density)
    {
        if (!std::isfinite(value))
        {
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc rVV10 received a nonfinite total density");
        }
    }
    for (double value: valence_density)
    {
        if (!std::isfinite(value))
        {
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc rVV10 received a nonfinite valence density");
        }
    }

    pw.real2recip(total_density.data(), density_g.data());
    XC_Functional::grad_rho(density_g.data(), gradient.data(), &pw, pw.tpiba);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        sigma[r] = gradient[r].x * gradient[r].x + gradient[r].y * gradient[r].y + gradient[r].z * gradient[r].z;
        field[r] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }

    std::vector<double> density = total_density;
    const double energy_ha = vdwxc_calculate(data_,
                                             density.data(),
                                             sigma.data(),
                                             dedn.data(),
                                             dedsigma.data());
    if (!std::isfinite(energy_ha))
    {
        ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc returned a nonfinite rVV10 energy");
    }

    for (int r = 0; r < pw.nrxx; ++r)
    {
        if (!std::isfinite(dedn[r]) || !std::isfinite(dedsigma[r]))
        {
            ModuleBase::WARNING_QUIT("LibvdwxcNonlocalXC", "libvdwxc returned a nonfinite rVV10 derivative");
        }
        const double coefficient = 2.0 * dedsigma[r];
        field[r].x = coefficient * gradient[r].x;
        field[r].y = coefficient * gradient[r].y;
        field[r].z = coefficient * gradient[r].z;
    }
    XC_Functional::grad_dot(field.data(), divergence.data(), &pw, pw.tpiba);

    const double energy_ry = 2.0 * energy_ha;
    std::vector<double> potential(pw.nrxx, 0.0);
    NonlocalXCResult result{energy_ry, 0.0, potential};
    const double dv = pw.omega / pw.nxyz;
    for (int r = 0; r < pw.nrxx; ++r)
    {
        // libvdwxc returns Hartree derivatives; ABACUS stores XC quantities in Ry.
        result.potential[r] = 2.0 * (dedn[r] - divergence[r]);
        if (!std::isfinite(divergence[r]) || !std::isfinite(result.potential[r]))
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

} // namespace vdw

#endif // __LIBVDWXC

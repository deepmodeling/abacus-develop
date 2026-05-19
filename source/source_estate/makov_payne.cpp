#include "makov_payne.h"

#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

namespace
{
constexpr double madelung_sc = 2.8373;
constexpr double madelung_fcc = 2.8883;
constexpr double madelung_bcc = 2.8885;

int cubic_ibrav(const std::string& latname)
{
    if (latname == "sc")
    {
        return 1;
    }
    if (latname == "fcc")
    {
        return 2;
    }
    if (latname == "bcc")
    {
        return 3;
    }
    return 0;
}

double madelung_constant(const int ibrav)
{
    if (ibrav == 1)
    {
        return madelung_sc;
    }
    if (ibrav == 2)
    {
        return madelung_fcc;
    }
    return madelung_bcc;
}

ModuleBase::Vector3<double> frac_to_cart_lat0(const UnitCell& ucell, const ModuleBase::Vector3<double>& frac)
{
    return frac.x * ucell.a1 + frac.y * ucell.a2 + frac.z * ucell.a3;
}

double wrap_delta(double x)
{
    x -= std::floor(x + 0.5);
    return x;
}

void add_average(const Charge& charge, const double* values, const int direction, const int length, std::vector<double>& ave)
{
    const ModulePW::PW_Basis* rhopw = charge.rhopw;
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        int index = 0;
        if (direction == 0)
        {
            index = ir / (rhopw->ny * rhopw->nplane);
        }
        else if (direction == 1)
        {
            const int ix = ir / (rhopw->ny * rhopw->nplane);
            index = ir / rhopw->nplane - ix * rhopw->ny;
        }
        else
        {
            index = ir % rhopw->nplane + rhopw->startz_current;
        }
        ave[index] += values[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ave.data(), length);
#endif
    const int surface = rhopw->nxyz / length;
    for (double& v : ave)
    {
        v /= static_cast<double>(surface);
    }
}

bool estimate_corrected_vacuum_level(const UnitCell& ucell,
                                      const Charge& charge,
                                      const double* v_elecstat,
                                      const ModuleBase::Vector3<double>& x0_frac,
                                      const double net_charge,
                                      double& vacuum_level_ev)
{
    if (v_elecstat == nullptr || charge.rhopw == nullptr || charge.rhopw->nxyz <= 0)
    {
        return false;
    }

    double vacuum[3] = {0.0, 0.0, 0.0};
    for (int dir = 0; dir < 3; ++dir)
    {
        std::vector<double> pos;
        pos.reserve(ucell.nat);
        for (int it = 0; it < ucell.ntype; ++it)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                pos.push_back(ucell.atoms[it].taud[ia][dir]);
            }
        }
        std::sort(pos.begin(), pos.end());
        for (std::size_t i = 1; i < pos.size(); ++i)
        {
            vacuum[dir] = std::max(vacuum[dir], pos[i] - pos[i - 1]);
        }
        if (!pos.empty())
        {
            vacuum[dir] = std::max(vacuum[dir], pos.front() + 1.0 - pos.back());
        }
    }

    int direction = 2;
    const double lengths[3] = {ucell.a1.norm() * ucell.lat0, ucell.a2.norm() * ucell.lat0, ucell.a3.norm() * ucell.lat0};
    vacuum[0] *= lengths[0];
    vacuum[1] *= lengths[1];
    vacuum[2] *= lengths[2];
    if (vacuum[0] > vacuum[direction])
    {
        direction = 0;
    }
    if (vacuum[1] > vacuum[direction])
    {
        direction = 1;
    }

    const ModulePW::PW_Basis* rhopw = charge.rhopw;
    int length = rhopw->nz;
    if (direction == 0)
    {
        length = rhopw->nx;
    }
    else if (direction == 1)
    {
        length = rhopw->ny;
    }

    std::vector<double> total_rho(rhopw->nrxx, 0.0);
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        total_rho[ir] = std::fabs(charge.rho[0][ir]);
    }
    if (PARAM.inp.nspin == 2)
    {
        for (int ir = 0; ir < rhopw->nrxx; ++ir)
        {
            total_rho[ir] += std::fabs(charge.rho[1][ir]);
        }
    }

    std::vector<double> rho_ave(length, 0.0);
    add_average(charge, total_rho.data(), direction, length, rho_ave);

    int min_index = 0;
    double min_value = 1.0e100;
    const double windows[7] = {0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};
    for (int i = 0; i < length; ++i)
    {
        double sum = 0.0;
        const int temp = i - 3 + length;
        for (int win = 0; win < 7; ++win)
        {
            sum += rho_ave[(temp + win) % length] * windows[win];
        }
        if (sum < min_value)
        {
            min_value = sum;
            min_index = i;
        }
    }

    std::vector<double> v_ave(length, 0.0);
    add_average(charge, v_elecstat, direction, length, v_ave);

    ModuleBase::Vector3<double> vac_frac = x0_frac;
    vac_frac[direction] = (static_cast<double>(min_index) + 0.5) / static_cast<double>(length);
    ModuleBase::Vector3<double> dfrac(wrap_delta(vac_frac.x - x0_frac.x),
                                      wrap_delta(vac_frac.y - x0_frac.y),
                                      wrap_delta(vac_frac.z - x0_frac.z));
    double r = frac_to_cart_lat0(ucell, dfrac).norm() * ucell.lat0;
    if (r < 1.0e-12)
    {
        r = 0.5 * std::min(lengths[0], std::min(lengths[1], lengths[2]));
    }

    vacuum_level_ev = (v_ave[min_index] + ModuleBase::e2 * net_charge / r) * ModuleBase::Ry_to_eV;
    return true;
}
} // namespace

namespace elecstate
{

MakovPayneResult makov_payne_correction(const UnitCell& ucell,
                                         const Charge& charge,
                                         const double* v_elecstat,
                                         bool print)
{
    const int ibrav = cubic_ibrav(PARAM.inp.latname);
    if (ibrav == 0)
    {
        ModuleBase::WARNING_QUIT("Makov-Payne", "Makov-Payne correction is defined only for cubic lattices: latname = sc, fcc, or bcc.");
    }
    if (charge.rhopw == nullptr || charge.rho == nullptr)
    {
        ModuleBase::WARNING_QUIT("Makov-Payne", "charge density is not available.");
    }

    double zion = 0.0;
    ModuleBase::Vector3<double> x0_frac(0.0, 0.0, 0.0);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const double zv = ucell.atoms[it].ncpp.zv;
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            zion += zv;
            x0_frac += ucell.atoms[it].taud[ia] * zv;
        }
    }
    if (zion <= 0.0)
    {
        ModuleBase::WARNING_QUIT("Makov-Payne",
                                 "total ionic valence must be positive to compute the Makov-Payne correction.");
    }
    x0_frac = x0_frac / zion;

    ModuleBase::Vector3<double> dipole_ion(0.0, 0.0, 0.0);
    ModuleBase::Vector3<double> quadrupole_ion(0.0, 0.0, 0.0);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const double zv = ucell.atoms[it].ncpp.zv;
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ModuleBase::Vector3<double> dfrac(wrap_delta(ucell.atoms[it].taud[ia].x - x0_frac.x),
                                              wrap_delta(ucell.atoms[it].taud[ia].y - x0_frac.y),
                                              wrap_delta(ucell.atoms[it].taud[ia].z - x0_frac.z));
            ModuleBase::Vector3<double> dr = frac_to_cart_lat0(ucell, dfrac) * ucell.lat0;
            dipole_ion += dr * zv;
            quadrupole_ion.x += zv * dr.x * dr.x;
            quadrupole_ion.y += zv * dr.y * dr.y;
            quadrupole_ion.z += zv * dr.z * dr.z;
        }
    }

    const ModulePW::PW_Basis* rhopw = charge.rhopw;
    const double dv = rhopw->omega / static_cast<double>(rhopw->nxyz);
    double electron_number = 0.0;
    ModuleBase::Vector3<double> dipole_el(0.0, 0.0, 0.0);
    ModuleBase::Vector3<double> quadrupole_el(0.0, 0.0, 0.0);

    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        const int ix = ir / (rhopw->ny * rhopw->nplane);
        const int iy = ir / rhopw->nplane - ix * rhopw->ny;
        const int iz = ir % rhopw->nplane + rhopw->startz_current;
        double rho = charge.rho[0][ir];
        if (PARAM.inp.nspin == 2)
        {
            rho += charge.rho[1][ir];
        }
        const double weight = rho * dv;
        electron_number += weight;
        ModuleBase::Vector3<double> frac((static_cast<double>(ix) + 0.5) / static_cast<double>(rhopw->nx),
                                         (static_cast<double>(iy) + 0.5) / static_cast<double>(rhopw->ny),
                                         (static_cast<double>(iz) + 0.5) / static_cast<double>(rhopw->nz));
        ModuleBase::Vector3<double> dfrac(wrap_delta(frac.x - x0_frac.x),
                                          wrap_delta(frac.y - x0_frac.y),
                                          wrap_delta(frac.z - x0_frac.z));
        ModuleBase::Vector3<double> dr = frac_to_cart_lat0(ucell, dfrac) * ucell.lat0;
        dipole_el += dr * weight;
        quadrupole_el.x += weight * dr.x * dr.x;
        quadrupole_el.y += weight * dr.y * dr.y;
        quadrupole_el.z += weight * dr.z * dr.z;
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(electron_number);
    Parallel_Reduce::reduce_pool(&dipole_el.x, 3);
    Parallel_Reduce::reduce_pool(&quadrupole_el.x, 3);
#endif

    MakovPayneResult result;
    result.charge = zion - electron_number;
    const ModuleBase::Vector3<double> dipole = dipole_ion - dipole_el;
    const ModuleBase::Vector3<double> quadrupole = quadrupole_ion - quadrupole_el;
    const double aa = quadrupole.x + quadrupole.y + quadrupole.z;
    const double bb = dipole.x * dipole.x + dipole.y * dipole.y + dipole.z * dipole.z;
    const double alpha = madelung_constant(ibrav);

    const double corr1 = -alpha / ucell.lat0 * result.charge * result.charge / 2.0 * ModuleBase::e2;
    const double corr2 = (2.0 / 3.0 * ModuleBase::PI) * (result.charge * aa - bb) / std::pow(ucell.lat0, 3)
                         * ModuleBase::e2;
    result.first_order = -corr1;
    result.second_order = -corr2;
    result.total = result.first_order + result.second_order;
    result.has_vacuum_level = estimate_corrected_vacuum_level(ucell,
                                                               charge,
                                                               v_elecstat,
                                                               x0_frac,
                                                               result.charge,
                                                               result.vacuum_level);

    if (print && GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running << std::setprecision(8) << std::fixed;
        GlobalV::ofs_running << "\n ********* MAKOV-PAYNE CORRECTION *********" << std::endl;
        GlobalV::ofs_running << " Makov-Payne correction with Madelung constant = " << std::setw(8) << alpha
                             << std::endl;
        GlobalV::ofs_running << " Makov-Payne correction " << std::setw(14) << result.first_order
                             << " Ry = " << std::setw(10) << result.first_order * ModuleBase::Ry_to_eV
                             << " eV (1st order, 1/a0)" << std::endl;
        GlobalV::ofs_running << "                         " << std::setw(14) << result.second_order
                             << " Ry = " << std::setw(10) << result.second_order * ModuleBase::Ry_to_eV
                             << " eV (2nd order, 1/a0^3)" << std::endl;
        GlobalV::ofs_running << "                         " << std::setw(14) << result.total
                             << " Ry = " << std::setw(10) << result.total * ModuleBase::Ry_to_eV
                             << " eV (total)" << std::endl;
        if (result.has_vacuum_level)
        {
            GlobalV::ofs_running << " Corrected vacuum level = " << std::setw(16) << result.vacuum_level << " eV"
                                 << std::endl;
        }
        else
        {
            GlobalV::ofs_running << " Corrected vacuum level was not evaluated because the electrostatic potential is unavailable."
                                 << std::endl;
        }
    }

    return result;
}

} // namespace elecstate

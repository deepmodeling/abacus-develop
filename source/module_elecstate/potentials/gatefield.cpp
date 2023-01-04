#include "gatefield.h"

#include "efield.h"
#include "module_base/timer.h"

namespace elecstate
{

double Gatefield::etotgatefield = 0.0;
double Gatefield::rho_surface;
double Gatefield::zgate = 0.5;
bool Gatefield::relax = false;
bool Gatefield::block = false;
double Gatefield::block_down = 0.45;
double Gatefield::block_up = 0.55;
double Gatefield::block_height = 0.1;

void Gatefield::add_gatefield(double *vltot,
                              const UnitCell &cell,
                              const ModulePW::PW_Basis *rho_basis,
                              const bool &linear,
                              const bool &quadratic)
{
    ModuleBase::TITLE("Gatefield", "add_gatefield");
    ModuleBase::timer::tick("Gatefield", "add_gatefield");

    //=======================================================
    // preparation for constants
    //=======================================================
    double latvec; // latvec along the efield direction
    double area; // surface area along the efield direction
    Efield::prepare(cell, latvec, area);

    double ion_charge = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        ion_charge += cell.atoms[it].na * cell.atoms[it].ncpp.zv;
    }
    rho_surface = -(GlobalV::nelec - ion_charge) / area * ModuleBase::TWO_PI;

    double block_size = block_up - block_down;

    //=======================================================
    // correction energy for gatefield
    //=======================================================
    double factor = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        double zval = cell.atoms[it].ncpp.zv;
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            double pos = cell.atoms[it].taud[ia][Efield::efield_dir];
            factor += zval * (mopopla(zgate, pos, true) + 1.0 / 6.0); // linear part
            factor += zval * mopopla(zgate, pos, false); // quadratic part
        }
    }
    etotgatefield
        = -ModuleBase::e2 * rho_surface * cell.lat0 / Efield::bmod * (factor + (GlobalV::nelec - ion_charge) / 12.0);

    GlobalV::ofs_running << "\n\n Adding charged plate to compensate the charge of the system" << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "prefactor of the potential (Ry a.u.)", rho_surface);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "position of the charged plate within cell", zgate);
    if (relax)
    {
        GlobalV::ofs_running << "Allow relaxation in specific direction" << std::endl;
    }
    if (block)
    {
        if (GlobalV::DIP_COR_FLAG)
        {
            GlobalV::ofs_running << "Adding potential to prevent charge spilling into region of the gate" << std::endl;
        }
        else
        {
            GlobalV::ofs_running << "Adding potential to prevent interaction between lower and upper part of unit cell"
                                 << std::endl;
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "block_size in units of unit cell length", block_size);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "block_height (Ry a.u.)", block_height);
    }

    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        int i = ir / (rho_basis->ny * rho_basis->nplane);
        int j = ir / rho_basis->nplane - i * rho_basis->ny;
        int k = ir % rho_basis->nplane + rho_basis->startz_current;
        double x = (double)i / rho_basis->nx;
        double y = (double)j / rho_basis->ny;
        double z = (double)k / rho_basis->nz;
        ModuleBase::Vector3<double> pos(x, y, z);
        double gatepos = pos[Efield::efield_dir];

        double value = 0;

        if (linear)
        {
            value += rho_surface * ModuleBase::e2 * (mopopla(zgate, gatepos, true) + 1.0 / 6.0) * cell.lat0
                     / Efield::bmod;
            if (block && gatepos >= block_down && gatepos <= block_up && !GlobalV::DIP_COR_FLAG)
            {
                if (gatepos - zgate <= -block_size / 2.0 * 0.9) // smooth increase within the first 10%
                {
                    value += block_height * (gatepos - zgate + block_size / 2.0) / (0.1 * block_size / 2.0);
                }
                else if (gatepos - zgate >= block_size / 2.0 * 0.9) // smooth decrease within the last 10%
                {
                    value += block_height * (gatepos - zgate - block_size / 2.0) / (-0.1 * block_size / 2.0);
                }
                else // block
                {
                    value += block_height;
                }
            }
            else if (block && gatepos >= block_down && gatepos <= block_up && GlobalV::DIP_COR_FLAG)
            {
                if (gatepos <= block_down + Efield::efield_pos_dec)
                {
                    value += (gatepos - block_down) / Efield::efield_pos_dec * block_height;
                }
                else if (gatepos >= block_up - Efield::efield_pos_dec)
                {
                    value += (block_up - gatepos) / Efield::efield_pos_dec * block_height;
                }
                else
                {
                    value += block_height;
                }
            }
        }

        if (quadratic)
        {
            value += rho_surface * ModuleBase::e2 * mopopla(zgate, gatepos, false) * cell.lat0 / Efield::bmod;
        }

        vltot[ir] += value;
    }

    ModuleBase::timer::tick("Gatefield", "add_gatefield");
}

double Gatefield::mopopla(double &zgate, double z, bool flag)
{
    while (z > 1.0)
        z -= 1.0;
    while (z < 0.0)
        z += 1.0;

    z -= zgate;

    if (z <= -0.5)
        z += 1.0;
    if (z >= 0.5)
        z -= 1.0;

    if (!flag)
    {
        return z * z;
    }
    else if (z <= 0)
    {
        return z;
    }
    else
    {
        return -z;
    }
}

void Gatefield::compute_force(const UnitCell &cell, ModuleBase::matrix &fgate)
{
    int iat = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        double zval = cell.atoms[it].ncpp.zv;
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            double pos = cell.atoms[it].taud[ia][Efield::efield_dir];
            while (pos > 1.0)
                pos -= 1.0;
            while (pos < 0.0)
                pos += 1.0;

            pos -= zgate;

            if (pos <= -0.5)
                pos += 1.0;
            if (pos >= 0.5)
                pos -= 1.0;

            double fac = 1;
            if (pos < 0)
                fac = -1;

            for (int jj = 0; jj < 3; ++jj)
            {
                fgate(iat, jj)
                    = -zval * ModuleBase::e2 * rho_surface * Efield::bvec[jj] / Efield::bmod * (fac - 2.0 * pos);
            }
            ++iat;
        }
    }
}

} // namespace elecstate

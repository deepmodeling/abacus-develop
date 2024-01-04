#include "efield.h"

#include "gatefield.h"
#include "module_base/constants.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"

namespace elecstate
{

double Efield::etotefield = 0.0;
double Efield::tot_dipole;
int Efield::efield_dir;
double Efield::efield_pos_max;
double Efield::efield_pos_dec;
double Efield::efield_amp;
double Efield::bvec[3];
double Efield::bmod;

Efield::Efield()
{
}

Efield::~Efield()
{
}

//=======================================================
// calculate dipole potential in surface calculations
//=======================================================
ModuleBase::matrix Efield::add_efield(const UnitCell& cell,
                                      const ModulePW::PW_Basis* rho_basis,
                                      const int& nspin,
                                      const double* const* const rho,
                                      surchem& solvent)
{
    ModuleBase::TITLE("Efield", "add_efield");
    ModuleBase::timer::tick("Efield", "add_efield");

    // set the parameters
    autoset(cell);

    double latvec; // latvec along the efield direction
    double area; // surface area along the efield direction
    prepare(cell, latvec, area);

    double ion_dipole = 0;
    double elec_dipole = 0;
    double induced_dipole = 0;

    if (GlobalV::DIP_COR_FLAG)
    {
        ion_dipole = cal_ion_dipole(cell, bmod);
        elec_dipole = cal_elec_dipole(cell, rho_basis, nspin, rho, bmod);
        tot_dipole = ion_dipole - elec_dipole;

        if (GlobalV::imp_sol)
        {
            induced_dipole = cal_induced_dipole(cell, rho_basis, solvent, bmod);
            tot_dipole += induced_dipole;
        }

        // energy correction
        etotefield = -ModuleBase::e2 * (efield_amp - 0.5 * tot_dipole) * tot_dipole * cell.omega / ModuleBase::FOUR_PI;
    }
    else
    {
        ion_dipole = cal_ion_dipole(cell, bmod);

        // energy correction
        etotefield = -ModuleBase::e2 * efield_amp * ion_dipole * cell.omega / ModuleBase::FOUR_PI;
    }

    const double length = (1.0 - efield_pos_dec) * latvec * cell.lat0;
    const double vamp = ModuleBase::e2 * (efield_amp - tot_dipole) * length;

    GlobalV::ofs_running << "\n\n Adding external electric field: " << std::endl;
    if (GlobalV::DIP_COR_FLAG)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Computed dipole along efield_dir", efield_dir);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Elec. dipole (Ry a.u.)", elec_dipole);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Ion dipole (Ry a.u.)", ion_dipole);
        if (GlobalV::imp_sol)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Induced dipole (Ry a.u.)", induced_dipole);
        }
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total dipole (Ry a.u.)", tot_dipole);
    }
    if (std::abs(efield_amp) > 0.0)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Amplitute of Efield (Hartree)", efield_amp);
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Potential amplitute (Ry)", vamp);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Total length (Bohr)", length);

    // dipole potential
    ModuleBase::matrix v(nspin, rho_basis->nrxx);
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        int i = ir / (rho_basis->ny * rho_basis->nplane);
        int j = ir / rho_basis->nplane - i * rho_basis->ny;
        int k = ir % rho_basis->nplane + rho_basis->startz_current;
        double x = (double)i / rho_basis->nx;
        double y = (double)j / rho_basis->ny;
        double z = (double)k / rho_basis->nz;
        ModuleBase::Vector3<double> pos(x, y, z);

        double saw = saw_function(efield_pos_max, efield_pos_dec, pos[efield_dir]);

        for (int is = 0; is < nspin0; is++)
        {
            v(is, ir) = saw;
        }
    }

    double fac = ModuleBase::e2 * (efield_amp - tot_dipole) * cell.lat0 / bmod;

    ModuleBase::timer::tick("Efield", "add_efield");
    return v * fac;
}

//=======================================================
// calculate dipole density in surface calculations
//=======================================================
double Efield::cal_ion_dipole(const UnitCell &cell, const double &bmod)
{
    double ion_dipole = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        double sum = 0;
        for (int ia = 0; ia < cell.atoms[it].na; ++ia)
        {
            sum += saw_function(efield_pos_max, efield_pos_dec, cell.atoms[it].taud[ia][efield_dir]);
        }
        ion_dipole += sum * cell.atoms[it].ncpp.zv;
    }

    if (GlobalV::GATE_FLAG && GlobalV::DIP_COR_FLAG)
    {
        double ion_charge = 0;
        for (int it = 0; it < cell.ntype; ++it)
        {
            ion_charge += cell.atoms[it].na * cell.atoms[it].ncpp.zv;
        }
        ion_dipole += (GlobalV::nelec - ion_charge) * saw_function(efield_pos_max, efield_pos_dec, Gatefield::zgate);
    }

    ion_dipole *= cell.lat0 / bmod * ModuleBase::FOUR_PI / cell.omega;

    return ion_dipole;
}

double Efield::cal_elec_dipole(const UnitCell& cell,
                               const ModulePW::PW_Basis* rho_basis,
                               const int& nspin,
                               const double* const* const rho,
                               const double& bmod)
{
    double elec_dipole = 0;
    const int nspin0 = (nspin == 2) ? 2 : 1;

    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        int i = ir / (rho_basis->ny * rho_basis->nplane);
        int j = ir / rho_basis->nplane - i * rho_basis->ny;
        int k = ir % rho_basis->nplane + rho_basis->startz_current;
        double x = (double)i / rho_basis->nx;
        double y = (double)j / rho_basis->ny;
        double z = (double)k / rho_basis->nz;
        ModuleBase::Vector3<double> pos(x, y, z);

        double saw = saw_function(efield_pos_max, efield_pos_dec, pos[efield_dir]);

        for (int is = 0; is < nspin0; is++)
        {
            elec_dipole += rho[is][ir] * saw;
        }
    }

    Parallel_Reduce::reduce_pool(elec_dipole);
    elec_dipole *= cell.lat0 / bmod * ModuleBase::FOUR_PI / rho_basis->nxyz;

    return elec_dipole;
}

double Efield::cal_induced_dipole(const UnitCell& cell,
                                  const ModulePW::PW_Basis* rho_basis,
                                  surchem& solvent,
                                  const double& bmod)
{
    double induced_dipole = 0;

    double *induced_rho = new double[rho_basis->nrxx];
    solvent.induced_charge(cell, rho_basis, induced_rho);

    for (int ir = 0; ir < rho_basis->nrxx; ++ir)
    {
        int i = ir / (rho_basis->ny * rho_basis->nplane);
        int j = ir / rho_basis->nplane - i * rho_basis->ny;
        int k = ir % rho_basis->nplane + rho_basis->startz_current;
        double x = (double)i / rho_basis->nx;
        double y = (double)j / rho_basis->ny;
        double z = (double)k / rho_basis->nz;
        ModuleBase::Vector3<double> pos(x, y, z);

        double saw = saw_function(efield_pos_max, efield_pos_dec, pos[efield_dir]);

        induced_dipole += induced_rho[ir] * saw;
    }

    Parallel_Reduce::reduce_pool(induced_dipole);
    induced_dipole *= cell.lat0 / bmod * ModuleBase::FOUR_PI / rho_basis->nxyz;

    return induced_dipole;
}

double Efield::saw_function(const double &a, const double &b, const double &x)
{
    assert(x >= 0);
    assert(x <= 1);

    const double fac = 1 - b;

    if (x <= a)
    {
        return x - a + 0.5 * fac;
    }
    else if (x > (a + b))
    {
        return x - a - 1 + 0.5 * fac;
    }
    else
    {
        return 0.5 * fac - fac * (x - a) / b;
    }
}

void Efield::compute_force(const UnitCell &cell, ModuleBase::matrix &fdip)
{
    if (GlobalV::DIP_COR_FLAG)
    {
        int iat = 0;
        for (int it = 0; it < cell.ntype; ++it)
        {
            for (int ia = 0; ia < cell.atoms[it].na; ++ia)
            {
                for (int jj = 0; jj < 3; ++jj)
                {
                    fdip(iat, jj)
                        = ModuleBase::e2 * (efield_amp - tot_dipole) * cell.atoms[it].ncpp.zv * bvec[jj] / bmod;
                }
                ++iat;
            }
        }
    }
    else
    {
        int iat = 0;
        for (int it = 0; it < cell.ntype; ++it)
        {
            for (int ia = 0; ia < cell.atoms[it].na; ++ia)
            {
                for (int jj = 0; jj < 3; ++jj)
                {
                    fdip(iat, jj) = ModuleBase::e2 * efield_amp * cell.atoms[it].ncpp.zv * bvec[jj] / bmod;
                }
                ++iat;
            }
        }
    }
}

void Efield::prepare(const UnitCell &cell, double &latvec, double &area)
{
    if (efield_dir == 0)
    {
        bvec[0] = cell.G.e11;
        bvec[1] = cell.G.e12;
        bvec[2] = cell.G.e13;
        latvec = cell.a1.norm();
        area = cross(cell.a2, cell.a3).norm() * cell.lat0 * cell.lat0;
    }
    else if (efield_dir == 1)
    {
        bvec[0] = cell.G.e21;
        bvec[1] = cell.G.e22;
        bvec[2] = cell.G.e23;
        latvec = cell.a2.norm();
        area = cross(cell.a3, cell.a1).norm() * cell.lat0 * cell.lat0;
    }
    else if (efield_dir == 2)
    {
        bvec[0] = cell.G.e31;
        bvec[1] = cell.G.e32;
        bvec[2] = cell.G.e33;
        latvec = cell.a3.norm();
        area = cross(cell.a1, cell.a2).norm() * cell.lat0 * cell.lat0;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Efield::prepare", "direction is wrong!");
    }
    bmod = sqrt(pow(bvec[0], 2) + pow(bvec[1], 2) + pow(bvec[2], 2));
}

void Efield::autoset(const UnitCell& cell)
{
    if (efield_pos_max == -1 || efield_pos_dec == -1)
    {
        // obtain the position of atoms along the efield direction
        std::vector<double> pos;
        for (int it = 0; it < cell.ntype; ++it)
        {
            for (int ia = 0; ia < cell.atoms[it].na; ++ia)
            {
                pos.push_back(cell.atoms[it].taud[ia][efield_dir]);
            }
        }

        // determine the vacuum region
        std::sort(pos.begin(), pos.end());
        double vacuum = 0.0;
        double center = 0.0;
        for (int i = 1; i < pos.size(); i++)
        {
            double diff = pos[i] - pos[i - 1];
            if (diff > vacuum)
            {
                vacuum = diff;
                center = (pos[i] + pos[i - 1]) / 2;
            }
        }

        // consider the periodic boundary condition
        double diff = pos[0] + 1 - pos[pos.size() - 1];
        if (diff > vacuum)
        {
            vacuum = diff;
            center = (pos[0] + pos[pos.size() - 1] + 1) / 2;
        }

        // set the parameters
        efield_pos_max = center - vacuum / 20;
        efield_pos_dec = vacuum / 10;
        while (efield_pos_max >= 1)
        {
            efield_pos_max -= 1;
        }
        while (efield_pos_max < 0)
        {
            efield_pos_max += 1;
        }

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Autoset efield_pos_max", efield_pos_max);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Autoset efield_pos_dec", efield_pos_dec);
    }
}

} // namespace elecstate
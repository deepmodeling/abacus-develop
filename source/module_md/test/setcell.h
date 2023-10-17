//==========================================================
// Used to genarate ucell
//==========================================================
#ifndef SETCELL_H
#define SETCELL_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"

Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

namespace ModuleIO
{
void print_force(std::ofstream& ofs_running,
                 const UnitCell& cell,
                 const std::string& name,
                 const ModuleBase::matrix& force,
                 bool ry = true)
{
}
void printstress_total(const ModuleBase::matrix& scs, bool ry)
{
}
} // namespace ModuleIO

class Setcell
{
  public:
    static void setupcell(UnitCell& ucell)
    {
        ucell.ntype = 1;

        ucell.atoms = new Atom[ucell.ntype];
        ucell.set_atom_flag = true;

        delete[] ucell.atom_label;
        delete[] ucell.atom_mass;
        ucell.atom_mass = new double[ucell.ntype];
        ucell.atom_label = new std::string[ucell.ntype];
        ucell.atom_mass[0] = 39.948;
        ucell.atom_label[0] = "Ar";

        ucell.lat0 = 1;
        ucell.lat0_angstrom = ucell.lat0 * 0.529177;
        ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
        ucell.tpiba2 = ucell.tpiba * ucell.tpiba;

        ucell.latvec.e11 = ucell.latvec.e22 = ucell.latvec.e33 = 10;
        ucell.latvec.e12 = ucell.latvec.e13 = ucell.latvec.e23 = 0;
        ucell.latvec.e21 = ucell.latvec.e31 = ucell.latvec.e32 = 0;

        ucell.a1.x = ucell.latvec.e11;
        ucell.a1.y = ucell.latvec.e12;
        ucell.a1.z = ucell.latvec.e13;

        ucell.a2.x = ucell.latvec.e21;
        ucell.a2.y = ucell.latvec.e22;
        ucell.a2.z = ucell.latvec.e23;

        ucell.a3.x = ucell.latvec.e31;
        ucell.a3.y = ucell.latvec.e32;
        ucell.a3.z = ucell.latvec.e33;

        ucell.nat = 4;
        ucell.atoms[0].na = 4;
        ucell.init_vel = 1;

        delete[] ucell.atoms[0].tau;
        delete[] ucell.atoms[0].taud;
        delete[] ucell.atoms[0].vel;
        delete[] ucell.atoms[0].mbl;
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[4];
        ucell.atoms[0].dis = new ModuleBase::Vector3<double>[4];
        ucell.atoms[0].taud = new ModuleBase::Vector3<double>[4];
        ucell.atoms[0].vel = new ModuleBase::Vector3<double>[4];
        ucell.atoms[0].mbl = new ModuleBase::Vector3<int>[4];
        ucell.atoms[0].mass = ucell.atom_mass[0];

        ucell.atoms[0].angle1 = new double[4];
        ucell.atoms[0].angle2 = new double[4];
        ucell.atoms[0].m_loc_ = new ModuleBase::Vector3<double>[4];

        ucell.atoms[0].taud[0].set(0.0, 0.0, 0.0);
        ucell.atoms[0].taud[1].set(0.52, 0.52, 0.0);
        ucell.atoms[0].taud[2].set(0.51, 0.0, 0.5);
        ucell.atoms[0].taud[3].set(0.0, 0.53, 0.5);
        ucell.atoms[0].vel[0].set(-0.000132080736364, 7.13429429835e-05, -1.40179977966e-05);
        ucell.atoms[0].vel[1].set(0.000153039878532, -0.000146533266608, 9.64491480698e-05);
        ucell.atoms[0].vel[2].set(-0.000133789480226, -3.0451038112e-06, -5.40998380137e-05);
        ucell.atoms[0].vel[3].set(0.000112830338059, 7.82354274358e-05, -2.83313122596e-05);
        for (int ia = 0; ia < 4; ++ia)
        {
            ucell.atoms[0].tau[ia] = ucell.atoms[0].taud[ia] * ucell.latvec;
            ucell.atoms[0].mbl[ia].set(1, 1, 1);
        }

        ucell.omega = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;

        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.GGT = ucell.G * ucell.GT;
        ucell.invGGT = ucell.GGT.Inverse();

        ucell.GT0 = ucell.latvec.Inverse();
        ucell.G0 = ucell.GT.Transpose();
        ucell.GGT0 = ucell.G * ucell.GT;
        ucell.invGGT0 = ucell.GGT.Inverse();

        ucell.set_iat2itia();
    };

    static void neighbor(Grid_Driver& grid_neigh, UnitCell& ucell)
    {
        GlobalV::SEARCH_RADIUS = 8.5 * ModuleBase::ANGSTROM_AU;
        INPUT.mdp.lj_rcut = 8.5 * ModuleBase::ANGSTROM_AU;
        INPUT.mdp.lj_epsilon = 0.01032 / ModuleBase::Hartree_to_eV;
        INPUT.mdp.lj_sigma = 3.405 * ModuleBase::ANGSTROM_AU;
        atom_arrange::search(1, GlobalV::ofs_running, grid_neigh, ucell, GlobalV::SEARCH_RADIUS, 0, 0);
    };

    static void parameters()
    {
        GlobalV::global_out_dir = "./";
        GlobalV::global_readin_dir = "./";
        GlobalV::SEARCH_RADIUS = 8.5 * ModuleBase::ANGSTROM_AU;
        GlobalV::CAL_STRESS = 1;

        INPUT.mdp.dump_virial = true;
        INPUT.mdp.dump_force = true;
        INPUT.mdp.dump_vel = true;
        INPUT.mdp.cal_stress = true;

        INPUT.mdp.md_restart = 0;
        INPUT.mdp.md_dt = 1;
        INPUT.mdp.md_tfirst = INPUT.mdp.md_tlast = 300;

        GlobalV::ESOLVER_TYPE = "lj";
        INPUT.mdp.lj_rcut = 8.5;
        INPUT.mdp.lj_epsilon = 0.01032;
        INPUT.mdp.lj_sigma = 3.405;

        INPUT.mdp.msst_direction = 2;
        INPUT.mdp.msst_qmass = 1;
        INPUT.mdp.msst_vel = 0;
        INPUT.mdp.msst_vis = 0;
        INPUT.mdp.msst_tscale = 0.01;

        INPUT.mdp.md_tfreq = 1;
        INPUT.mdp.md_tchain = 4;
        INPUT.mdp.md_pfreq = 1;
        INPUT.mdp.md_pchain = 4;

        INPUT.mdp.md_damp = 1;

        INPUT.mdp.md_nraise = 2;
        INPUT.mdp.md_tolerance = 0;
    };
};

#endif
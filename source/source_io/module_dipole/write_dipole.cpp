#include "source_base/parallel_reduce.h"
#include "source_io/module_dipole/dipole_io.h"

// Calculate and write dipole moment data for RT-TDDFT calculations
// 
// Dipole moment is a measure of the separation of positive and negative charges
// in a system. In electronic structure calculations, we compute three components:
// 
// 1. Electronic dipole moment (P_elec):
//    Formula: P_elec = -int(r * rho(r) dr)
//    where rho(r) is the electron density
//
// 2. Ionic dipole moment (P_ion):
//    Formula: P_ion = sum_atoms (Z_v * r_atom)
//    where Z_v is the valence charge of the atom
//
// 3. Total dipole moment (P_tot):
//    Formula: P_tot = P_elec + P_ion
//
// The total dipole moment norm is |P_tot| = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
void ModuleIO::write_dipole(const UnitCell& ucell,
                            const double* rho_save,
                            const ModulePW::PW_Basis* rhopw,
                            const int& is,
                            const int& istep,
                            const std::string& fn,
                            std::ofstream& ofs_running,
                            const int& precision)
{
    ModuleBase::TITLE("ModuleIO", "write_dipole");

    time_t start, end;
    std::ofstream ofs;

    if (GlobalV::MY_RANK == 0)
    {
        start = time(NULL);

        ofs.open(fn.c_str(), std::ofstream::app);
        if (!ofs)
        {
            ModuleBase::WARNING("ModuleIO", "Can't create Dipole File!");
        }
    }

    ofs_running << " Write dipole data to file: " << fn << std::endl;

    // Calculate modulus of reciprocal lattice vectors
    // bmod[i] = |b_i| where b_i are reciprocal lattice vectors
    // Used for coordinate transformation from fractional to Cartesian
    double bmod[3];
    for (int i = 0; i < 3; i++)
    {
        bmod[i] = prepare(ucell, i);
        if (bmod[i] < 1e-10)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", "bmod[" + std::to_string(i) + "] is zero or too small!");
        }
    }

    // Validate grid dimensions to prevent division by zero
    if (rhopw->nx == 0 || rhopw->ny == 0 || rhopw->nz == 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", "Grid dimensions nx, ny, or nz cannot be zero!");
    }

    if (rhopw->nxyz == 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", "Total grid points nxyz cannot be zero!");
    }

    if (rhopw->nplane == 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", "nplane cannot be zero!");
    }

    // Calculate electronic dipole moment
    // P_elec = -int(r * rho(r) dr)
    // Discretized: P_elec[i] = -sum(r_grid[i] * rho(r_grid)) * (omega / nxyz)
    double dipole_elec[3] = {0.0, 0.0, 0.0};

    // Loop over local grid points (parallel decomposition)
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        // Convert 1D index to 3D indices
        int i = ir / (rhopw->ny * rhopw->nplane);
        int j = ir / rhopw->nplane - i * rhopw->ny;
        int k = ir % rhopw->nplane + rhopw->startz_current;
        
        // Convert to fractional coordinates: r_i = i / N_i
        double x = (double)i / rhopw->nx;
        double y = (double)j / rhopw->ny;
        double z = (double)k / rhopw->nz;

        // Accumulate: P_elec -= rho * r (negative sign from electron charge)
        dipole_elec[0] -= rho_save[ir] * x;
        dipole_elec[1] -= rho_save[ir] * y;
        dipole_elec[2] -= rho_save[ir] * z;
    }

    // Reduce across MPI processes to get global sum
    Parallel_Reduce::reduce_pool(dipole_elec[0]);
    Parallel_Reduce::reduce_pool(dipole_elec[1]);
    Parallel_Reduce::reduce_pool(dipole_elec[2]);

    // Convert to Cartesian coordinates and normalize
    // Conversion factor: lat0 / bmod[i] transforms fractional to Cartesian
    // Volume normalization: omega / nxyz accounts for grid spacing
    for (int i = 0; i < 3; ++i)
    {
        dipole_elec[i] *= ucell.lat0 / bmod[i] * ucell.omega / rhopw->nxyz;
    }

    // Output electronic dipole moment
    ofs_running << " Electronic dipole moment" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_elec_x(t)", dipole_elec[0]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_elec_y(t)", dipole_elec[1]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_elec_z(t)", dipole_elec[2]);

    // Write to file: step index and three dipole components
    ofs << std::setprecision(precision) << istep+1 
        << " " << dipole_elec[0] 
	    << " " << dipole_elec[1] 
	    << " " << dipole_elec[2] << std::endl;

    // Calculate ionic dipole moment
    // P_ion = sum_{atom_types} sum_{atoms} (Z_v * tau)
    // where tau is the atomic position in fractional coordinates
    double dipole_ion[3] = {0.0};
    double dipole_sum = 0.0;

    for (int i = 0; i < 3; ++i)
    {
        for (int it = 0; it < ucell.ntype; ++it)
        {
            double sum = 0;
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                sum += ucell.atoms[it].taud[ia][i];
            }
            dipole_ion[i] += sum * ucell.atoms[it].ncpp.zv;
        }
        // Convert to Cartesian coordinates
        dipole_ion[i] *= ucell.lat0 / bmod[i];
    }

    // Output ionic dipole moment
    ofs_running << " Ionic dipole moment" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_ion_x(t)", dipole_ion[0]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_ion_y(t)", dipole_ion[1]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_ion_z(t)", dipole_ion[2]);

    // Calculate total dipole moment
    // P_tot = P_elec + P_ion
    double dipole[3] = {0.0};
    for (int i = 0; i < 3; ++i)
    {
        dipole[i] = dipole_ion[i] + dipole_elec[i];
    }

    // Output total dipole moment
    ofs_running << " Total dipole moment" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_tot_x(t)", dipole[0]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_tot_y(t)", dipole[1]);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_tot_z(t)", dipole[2]);

    // Calculate and output total dipole moment norm
    // |P_tot| = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
    dipole_sum = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);

    ofs_running << " Total dipole moment norm" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "|P_tot(t)|", dipole_sum);

    if (GlobalV::MY_RANK == 0)
    {
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_dipole", start, end);
        ofs.close();
    }

    return;
}

// Calculate the modulus of a reciprocal lattice vector
// Input: cell - unit cell containing reciprocal lattice G
//        dir - direction index (0=x, 1=y, 2=z)
// Output: bmod - |b_dir| where b_dir is the reciprocal lattice vector
double ModuleIO::prepare(const UnitCell& cell, int& dir)
{
    double bvec[3] = {0.0};
    double bmod = 0.0;
    if (dir == 0)
    {
        bvec[0] = cell.G.e11;
        bvec[1] = cell.G.e12;
        bvec[2] = cell.G.e13;
    }
    else if (dir == 1)
    {
        bvec[0] = cell.G.e21;
        bvec[1] = cell.G.e22;
        bvec[2] = cell.G.e23;
    }
    else if (dir == 2)
    {
        bvec[0] = cell.G.e31;
        bvec[1] = cell.G.e32;
        bvec[2] = cell.G.e33;
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::prepare", "direction is wrong!");
    }
    bmod = sqrt(pow(bvec[0], 2) + pow(bvec[1], 2) + pow(bvec[2], 2));
    return bmod;
}
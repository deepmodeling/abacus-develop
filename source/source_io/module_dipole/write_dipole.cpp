#include "source_base/parallel_reduce.h"
#include "source_io/module_dipole/dipole_io.h"

namespace ModuleIO
{

// Helper function to print dipole moment components
// name: descriptive name for the dipole moment type
// px, py, pz: dipole moment components in x, y, z directions
void print_dipole_moment(std::ofstream& ofs_running, const std::string& name,
                       double px, double py, double pz)
{
    ofs_running << " " << name << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_x(t)", px);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_y(t)", py);
    ModuleBase::GlobalFunc::OUT(ofs_running, "P_z(t)", pz);
}

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
//    Formula: P_ion = sum_{atom_types} sum_{atoms} (Z_v * tau)
//    where Z_v is the valence charge and tau is atomic position
//
// 3. Total dipole moment (P_tot):
//    Formula: P_tot = P_elec + P_ion
//
// The total dipole moment norm is |P_tot| = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
//
// Parameters:
// - ucell: unit cell containing atomic structure and lattice information
// - rho_save: electron density on real-space grid
// - rhopw: plane wave basis information including grid dimensions
// - istep: current time step
// - fn: output file name
// - ofs_running: output stream for logging
// - precision: floating-point precision for output
void write_dipole(const UnitCell& ucell,
                  const double* rho_save,
                  const ModulePW::PW_Basis* rhopw,
                  const int& istep,
                  const std::string& fn,
                  std::ofstream& ofs_running,
                  const int& precision)
{
    ModuleBase::TITLE("ModuleIO", "write_dipole");

    time_t start, end;
    std::ofstream ofs;

    // Open output file on master process only
    if (GlobalV::MY_RANK == 0)
    {
        start = time(NULL);
        ofs.open(fn.c_str(), std::ofstream::app);
        if (!ofs)
        {
            ModuleBase::WARNING_QUIT("ModuleIO", "Cannot create dipole file: " + fn);
        }
    }

    ofs_running << " Write dipole data to file: " << fn << std::endl;

    // Calculate modulus of reciprocal lattice vectors
    // bmod[i] = |b_i| where b_i are reciprocal lattice vectors
    // Used for coordinate transformation from fractional to Cartesian
    double small_value = 1e-10;
    double bmod[3];
    
    for (int i = 0; i < 3; ++i)
    {
        bmod[i] = prepare(ucell, i);
        if (bmod[i] < small_value)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", 
                "bmod[" + std::to_string(i) + "] is too small or zero");
        }
    }

    // Validate grid dimensions to prevent division by zero
    if (rhopw->nx == 0 || rhopw->ny == 0 || rhopw->nz == 0 || 
        rhopw->nxyz == 0 || rhopw->nplane == 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", 
            "Invalid grid parameters: nx, ny, nz, nxyz, or nplane is zero");
    }

    // Calculate electronic dipole moment
    // P_elec = -int(r * rho(r) dr)
    // Discretized: P_elec[i] = -sum(r_grid[i] * rho(r_grid)) * (omega / nxyz)
    double dipole_elec[3] = {0.0, 0.0, 0.0};
    
    // Precompute inverse grid dimensions for performance
    double inv_nx = 1.0 / static_cast<double>(rhopw->nx);
    double inv_ny = 1.0 / static_cast<double>(rhopw->ny);
    double inv_nz = 1.0 / static_cast<double>(rhopw->nz);

    // Loop over local grid points (parallel decomposition with OpenMP)
    // Use reduction for thread-safe accumulation
    #pragma omp parallel for reduction(-:dipole_elec[:3]) schedule(static)
    for (int ir = 0; ir < rhopw->nrxx; ++ir)
    {
        // Convert 1D index to 3D indices
        int i = ir / (rhopw->ny * rhopw->nplane);
        int j = ir / rhopw->nplane - i * rhopw->ny;
        int k = ir % rhopw->nplane + rhopw->startz_current;
        
        // Convert to fractional coordinates: r_i = i / N_i
        // Using multiplication instead of division for better performance
        double x = static_cast<double>(i) * inv_nx;
        double y = static_cast<double>(j) * inv_ny;
        double z = static_cast<double>(k) * inv_nz;

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
    double vol_factor = ucell.omega / static_cast<double>(rhopw->nxyz);
    for (int i = 0; i < 3; ++i)
    {
        dipole_elec[i] *= ucell.lat0 / bmod[i] * vol_factor;
    }

    // Output electronic dipole moment
    print_dipole_moment(ofs_running, "Electronic dipole moment",
                      dipole_elec[0], dipole_elec[1], dipole_elec[2]);

    // Write to file: step index and three dipole components
    ofs << std::setprecision(precision) << istep + 1 
        << " " << dipole_elec[0] 
        << " " << dipole_elec[1] 
        << " " << dipole_elec[2] << std::endl;

    // Calculate ionic dipole moment
    // P_ion = sum_{atom_types} sum_{atoms} (Z_v * tau)
    // where tau is the atomic position in fractional coordinates
    double dipole_ion[3] = {0.0};
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
    print_dipole_moment(ofs_running, "Ionic dipole moment",
                      dipole_ion[0], dipole_ion[1], dipole_ion[2]);

    // Calculate total dipole moment
    // P_tot = P_elec + P_ion
    double dipole[3] = {0.0};
    for (int i = 0; i < 3; ++i)
    {
        dipole[i] = dipole_ion[i] + dipole_elec[i];
    }

    // Output total dipole moment
    print_dipole_moment(ofs_running, "Total dipole moment",
                      dipole[0], dipole[1], dipole[2]);

    // Calculate and output total dipole moment norm
    // |P_tot| = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
    double dipole_sum = sqrt(dipole[0] * dipole[0] + 
                             dipole[1] * dipole[1] + 
                             dipole[2] * dipole[2]);
    ofs_running << " Total dipole moment norm" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "|P_tot(t)|", dipole_sum);

    // Close file and report timing on master process
    if (GlobalV::MY_RANK == 0)
    {
        end = time(NULL);
        ModuleBase::GlobalFunc::OUT_TIME("write_dipole", start, end);
        ofs.close();
    }
}

// Calculate the modulus of a reciprocal lattice vector
// Input: 
//   cell - unit cell containing reciprocal lattice G
//   dir - direction index (0=x, 1=y, 2=z)
// Output: 
//   bmod - |b_dir| where b_dir is the reciprocal lattice vector
double prepare(const UnitCell& cell, int& dir)
{
    double bvec[3] = {0.0};
    
    // Select the appropriate reciprocal lattice vector components
    switch (dir)
    {
        case 0:  // x-direction
            bvec[0] = cell.G.e11;
            bvec[1] = cell.G.e12;
            bvec[2] = cell.G.e13;
            break;
        case 1:  // y-direction
            bvec[0] = cell.G.e21;
            bvec[1] = cell.G.e22;
            bvec[2] = cell.G.e23;
            break;
        case 2:  // z-direction
            bvec[0] = cell.G.e31;
            bvec[1] = cell.G.e32;
            bvec[2] = cell.G.e33;
            break;
        default:
            ModuleBase::WARNING_QUIT("ModuleIO::prepare", "Invalid direction index");
    }
    
    // Calculate and return the magnitude of the reciprocal lattice vector
    return sqrt(bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
}

} // namespace ModuleIO
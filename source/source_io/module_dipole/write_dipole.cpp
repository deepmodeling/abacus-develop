#include "source_base/parallel_reduce.h"
#include "source_io/module_dipole/dipole_io.h"

namespace ModuleIO
{

// Helper function to print dipole moment components
// name: descriptive name for the dipole moment type
// px, py, pz: dipole moment components in x, y, z directions
void print_dipole_moment(std::ofstream& ofs_running, const std::string& name, double px, double py, double pz)
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
// - rho: current electron density on the real-space grid
// - rhopw: plane wave basis information including grid dimensions
// - istep: current time step
// - fn: output file name
// - ofs_running: output stream for logging
// - precision: floating-point precision for output
void write_dipole(const UnitCell& ucell,
                  const double* rho,
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

    // Validate grid dimensions to prevent division by zero
    if (rhopw->nx == 0 || rhopw->ny == 0 || rhopw->nz == 0 || rhopw->nxyz == 0 || rhopw->nplane == 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dipole", "Invalid grid parameters: nx, ny, nz, nxyz, or nplane is zero");
    }

    // Accumulate the electronic first moment in fractional coordinates.
    double dipole_elec_frac[3] = {0.0, 0.0, 0.0};

    // Precompute inverse grid dimensions for performance
    double inv_nx = 1.0 / static_cast<double>(rhopw->nx);
    double inv_ny = 1.0 / static_cast<double>(rhopw->ny);
    double inv_nz = 1.0 / static_cast<double>(rhopw->nz);

// Loop over local grid points (parallel decomposition with OpenMP)
// Use reduction for thread-safe accumulation
#pragma omp parallel for reduction(- : dipole_elec_frac[ : 3]) schedule(static)
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

        // The negative sign accounts for the electron charge.
        dipole_elec_frac[0] -= rho[ir] * x;
        dipole_elec_frac[1] -= rho[ir] * y;
        dipole_elec_frac[2] -= rho[ir] * z;
    }

    // Reduce across MPI processes to get global sum
    Parallel_Reduce::reduce_pool(dipole_elec_frac[0]);
    Parallel_Reduce::reduce_pool(dipole_elec_frac[1]);
    Parallel_Reduce::reduce_pool(dipole_elec_frac[2]);

    // Lattice vectors are rows: r_cart = r_frac * latvec * lat0.
    // Apply the full matrix to retain shear, rotation, and axis signs.
    const double vol_factor = ucell.omega / static_cast<double>(rhopw->nxyz);
    const ModuleBase::Vector3<double> elec_frac(dipole_elec_frac[0], dipole_elec_frac[1], dipole_elec_frac[2]);
    const ModuleBase::Vector3<double> dipole_elec = elec_frac * ucell.latvec * (ucell.lat0 * vol_factor);

    // Output electronic dipole moment
    print_dipole_moment(ofs_running, "Electronic dipole moment", dipole_elec[0], dipole_elec[1], dipole_elec[2]);

    // Write to file: step index and three dipole components
    ofs << std::setprecision(precision) << istep + 1 << " " << dipole_elec[0] << " " << dipole_elec[1] << " " << dipole_elec[2]
        << std::endl;

    // Calculate ionic dipole moment
    // Accumulate Z_v * taud, then transform to Cartesian coordinates.
    ModuleBase::Vector3<double> ion_frac(0.0, 0.0, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        for (int it = 0; it < ucell.ntype; ++it)
        {
            double sum = 0;
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                sum += ucell.atoms[it].taud[ia][i];
            }
            ion_frac[i] += sum * ucell.atoms[it].ncpp.zv;
        }
    }
    const ModuleBase::Vector3<double> dipole_ion = ion_frac * ucell.latvec * ucell.lat0;

    // Output ionic dipole moment
    print_dipole_moment(ofs_running, "Ionic dipole moment", dipole_ion[0], dipole_ion[1], dipole_ion[2]);

    // Calculate total dipole moment
    // P_tot = P_elec + P_ion
    double dipole[3] = {0.0};
    for (int i = 0; i < 3; ++i)
    {
        dipole[i] = dipole_ion[i] + dipole_elec[i];
    }

    // Output total dipole moment
    print_dipole_moment(ofs_running, "Total dipole moment", dipole[0], dipole[1], dipole[2]);

    // Calculate and output total dipole moment norm
    // |P_tot| = sqrt(P_tot_x^2 + P_tot_y^2 + P_tot_z^2)
    double dipole_sum = sqrt(dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2]);
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

} // namespace ModuleIO

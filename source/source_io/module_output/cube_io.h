#ifndef CUBE_IO_H
#define CUBE_IO_H
#include "source_cell/unitcell.h"
#include "spin_tag.h"

#include <fstream>
#include <string>
#include <vector>
class Parallel_Grid;

namespace ModuleIO
{
/// Build a descriptive label for the data channel at spin index is.
/// For nspin==2 and is>=0: appends " (spin up  )" or " (spin down)".
/// For nspin==4 and is>0: returns mag_base with the axis letter appended,
///     e.g. mag_base "magnetization density m" yields "magnetization density mx".
/// For nspin==1, nspin==4 with is==0, or spin-summed calls (is<0): no suffix.
///
/// @param base label of the scalar channel (is == 0), e.g. "charge density".
/// @param mag_base label of the nspin==4 magnetic channels (is = 1..3),
///     ending with the component letter: "m" for magnetization densities
///     (rho[1..3]), "b" for magnetic fields (the B_xc components of the
///     effective potential).
inline std::string make_data_desc(const std::string& base, const std::string& mag_base, const int is, const int nspin)
{
    std::string desc = base + make_spin_tag(is, nspin);
    if (nspin == 4 && is > 0)
    {
        desc = mag_base;
        if (is == 1)
        {
            desc += "x";
        }
        else if (is == 2)
        {
            desc += "y";
        }
        else if (is == 3)
        {
            desc += "z";
        }
    }
    return desc;
}

/// read volumetric data from .cube file into the parallel distributed grid.
bool read_vdata_palgrid(const Parallel_Grid& pgrid,
                        const int my_rank,
                        std::ostream& ofs_running,
                        const std::string& fn,
                        double* const data,
                        const int nat);

/// write volumetric data on the parallized grid into a .cube file
void write_vdata_palgrid(const Parallel_Grid& pgrid,
                         const double* const data,
                         const int is,
                         const int nspin,
                         const int iter,
                         const std::string& fn,
                         const double ef,
                         const UnitCell* const ucell,
                         const int precision,
                         const int out_fermi,
                         const bool two_fermi,
                         const bool reduce_all_pool,
                         std::ofstream& ofs_running,
                         const std::string& data_desc);

/// read the full data from a cube file
bool read_cube(const std::string& file,
               std::vector<std::string>& comment,
               int& natom,
               std::vector<double>& origin,
               int& nx,
               int& ny,
               int& nz,
               std::vector<double>& dx,
               std::vector<double>& dy,
               std::vector<double>& dz,
               std::vector<int>& atom_type,
               std::vector<double>& atom_charge,
               std::vector<std::vector<double>>& atom_pos,
               std::vector<double>& data);

/// write a cube file
void write_cube(const std::string& file,
                const std::vector<std::string>& comment,
                const int& natom,
                const std::vector<double>& origin,
                const int& nx,
                const int& ny,
                const int& nz,
                const std::vector<double>& dx,
                const std::vector<double>& dy,
                const std::vector<double>& dz,
                const std::vector<int>& atom_type,
                const std::vector<double>& atom_charge,
                const std::vector<std::vector<double>>& atom_pos,
                const std::vector<double>& data,
                const int precision,
                const int ndata_line = 6);

/**
 * @brief The trilinear interpolation method
 *
 * Trilinear interpolation is a method for interpolating grid data in 3D space.
 * It estimates the value at a given position by interpolating the data along the three adjacent points.
 *
 * Specifically, for 3D grid data, trilinear interpolation requires determining the eight data points that are
 * closest to the point where the estimation is required. These data points form a cube, with vertices at
 * (x0,y0,z0), (x0+1,y0,z0), (x0,y0+1,z0), (x0+1,y0+1,z0), (x0,y0,z0+1), (x0+1,y0,z0+1), (x0,y0+1,z0+1) and
 * (x0+1,y0+1,z0+1). Here, (x0,y0,z0) is the data point closest to the estimation point and has coordinate
 * values less than those of the estimation point.
 *
 * For the estimation location (x,y,z), its estimated value in the grid can be calculated using the following
 * formula: f(x,y,z) = f000(1-dx)(1-dy)(1-dz) + f100dx(1-dy)(1-dz) + f010(1-dx)dy(1-dz) +
 *             f001(1-dx)(1-dy)dz + f101dx(1-dy)dz + f011(1-dx)dydz +
 *             f110dxdy(1-dz) + f111dxdydz
 * where fijk represents the data value at vertex i,j,k of the cube, and dx = x - x0, dy = y - y0, dz = z - z0
 * represent the distance between the estimation point and the closest data point in each of the three
 * directions, divided by the grid spacing. Here, it is assumed that the grid spacing is equal and can be
 * omitted during computation.
 *
 * @param data_in the input data of size nxyz_read
 * @param nx_read nx read from file
 * @param ny_read ny read from file
 * @param nz_read nz read from file
 * @param nx the dimension of grids along x
 * @param ny the dimension of grids along y
 * @param nz the dimension of grids along z
 * @param data_out the interpolated results of size nxyz
 */
void trilinear_interpolate(const double* const data_in,
                           const int& nx_read,
                           const int& ny_read,
                           const int& nz_read,
                           const int& nx,
                           const int& ny,
                           const int& nz,
                           double* data_out);
} // namespace ModuleIO

#endif

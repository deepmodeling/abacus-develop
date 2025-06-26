#ifndef CAL_LDOS_H
#define CAL_LDOS_H

#include "source_estate/elecstate_lcao.h"
#include "source_estate/elecstate_pw.h"

namespace ModuleIO
{
template <typename T>
class Cal_ldos
{
  public:
    Cal_ldos(){};
    ~Cal_ldos(){};

    static void cal_ldos_lcao(const elecstate::ElecStateLCAO<T>* pelec,
                              const psi::Psi<T>& psi,
                              const Parallel_Grid& pgrid,
                              const UnitCell& ucell);
}; // namespace Cal_ldos

void cal_ldos_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                 const psi::Psi<std::complex<double>>& psi,
                 const Parallel_Grid& pgrid,
                 const UnitCell& ucell);

void stm_mode_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                 const psi::Psi<std::complex<double>>& psi,
                 const Parallel_Grid& pgrid,
                 const UnitCell& ucell);

void ldos_mode_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                  const psi::Psi<std::complex<double>>& psi,
                  const Parallel_Grid& pgrid,
                  const UnitCell& ucell);

/*
 * @brief Get grid points and shifts for interpolation.
 *
 * @param start The starting point of the line.
 * @param end The ending point of the line.
 * @param npoints The number of points in the line.
 * @param nx The dimension of 3D grids in the x direction.
 * @param ny The dimension of 3D grids in the y direction.
 * @param nz The dimension of 3D grids in the z direction.
 * @param points The grid index that the points in the line are placed in.
 * @param shifts The shifts along three directions due to the points are not on the grid exactly.
 */
void get_grid_points(const std::vector<double>& start,
                     const std::vector<double>& end,
                     const int& npoints,
                     const int& nx,
                     const int& ny,
                     const int& nz,
                     std::vector<std::vector<int>>& points,
                     std::vector<std::vector<double>>& shifts);

/*
 * @brief Trilinear interpolation of data on a 3D grid.
 *
 * @param points The grid points for interpolation.
 * @param shifts The shifts for interpolation.
 * @param pgrid The parallel grid object.
 * @param data The data to be interpolated.
 * @param results The results of the interpolation.
 */
void trilinear_interpolate(const std::vector<std::vector<int>>& points,
                           const std::vector<std::vector<double>>& shifts,
                           const Parallel_Grid& pgrid,
                           const std::vector<double>& data,
                           std::vector<double>& results);

} // namespace ModuleIO

#endif // CAL_LDOS_H
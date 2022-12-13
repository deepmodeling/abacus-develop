#ifndef MATHZONE_H
#define MATHZONE_H

#include "global_function.h"
#include "matrix3.h"
#include "vector3.h"
#include "realarray.h"

#include <cassert>
#include <complex>
#include <map>
#include <vector>
#include <array>
namespace ModuleBase
{

/**
 * @brief atomic coordinates conversion functions
 *
 */
class Mathzone
{
  public:
    Mathzone();
    ~Mathzone();

  public:
    /**
     * @brief Pointwise product of two vectors with same size
     *
     * @tparam Type
     * @param[in] f1
     * @param[in] f2
     * @return std::vector<Type>
     * @author Peize Lin (2016-08-03)
     */
    template <typename Type>
    static std::vector<Type> Pointwise_Product(const std::vector<Type> &f1, const std::vector<Type> &f2)
    {
        assert(f1.size() == f2.size());
        std::vector<Type> f(f1.size());
        for (int ir = 0; ir != f.size(); ++ir)
            f[ir] = f1[ir] * f2[ir];
        return f;
    }

    /**
     * @brief change direct coordinate (dx,dy,dz) to
     * Cartesian coordinate (cx,cy,cz), (dx,dy,dz) = (cx,cy,cz) * R
     *
     * @param[in] dx Direct coordinats
     * @param[in] dy
     * @param[in] dz
     * @param[in] R11 Lattice vector matrix R_ij: i_row, j_column
     * @param[in] R12
     * @param[in] R13
     * @param[in] R21
     * @param[in] R22
     * @param[in] R23
     * @param[in] R31
     * @param[in] R32
     * @param[in] R33
     * @param[out] cx Cartesian coordinats
     * @param[out] cy
     * @param[out] cz 
     */
    static inline void Direct_to_Cartesian(const double &dx,
                                           const double &dy,
                                           const double &dz,
                                           const double &R11,
                                           const double &R12,
                                           const double &R13,
                                           const double &R21,
                                           const double &R22,
                                           const double &R23,
                                           const double &R31,
                                           const double &R32,
                                           const double &R33,
                                           double &cx,
                                           double &cy,
                                           double &cz)
    {
        ModuleBase::Matrix3 lattice_vector;
        ModuleBase::Vector3<double> direct_vec, cartesian_vec;
        lattice_vector.e11 = R11;
        lattice_vector.e12 = R12;
        lattice_vector.e13 = R13;
        lattice_vector.e21 = R21;
        lattice_vector.e22 = R22;
        lattice_vector.e23 = R23;
        lattice_vector.e31 = R31;
        lattice_vector.e32 = R32;
        lattice_vector.e33 = R33;

        direct_vec.x = dx;
        direct_vec.y = dy;
        direct_vec.z = dz;

        cartesian_vec = direct_vec * lattice_vector;
        cx = cartesian_vec.x;
        cy = cartesian_vec.y;
        cz = cartesian_vec.z;
        return;
    }

    /**
     * @brief Change Cartesian coordinate (cx,cy,cz) to
     * direct coordinate (dx,dy,dz), (cx,cy,cz) = (dx,dy,dz) * R^(-1)
     *
     * @param[in] cx Cartesian coordinats
     * @param[in] cy
     * @param[in] cz
     * @param[in] R11 Lattice vector matrix R_ij: i_row, j_column
     * @param[in] R12
     * @param[in] R13
     * @param[in] R21
     * @param[in] R22
     * @param[in] R23
     * @param[in] R31
     * @param[in] R32
     * @param[in] R33
     * @param[out] dx Direct coordinats
     * @param[out] dy
     * @param[out] dz
     */
    static inline void Cartesian_to_Direct(const double &cx,
                                           const double &cy,
                                           const double &cz,
                                           const double &R11,
                                           const double &R12,
                                           const double &R13,
                                           const double &R21,
                                           const double &R22,
                                           const double &R23,
                                           const double &R31,
                                           const double &R32,
                                           const double &R33,
                                           double &dx,
                                           double &dy,
                                           double &dz)
    {
        ModuleBase::Matrix3 lattice_vector, inv_lat;
        lattice_vector.e11 = R11;
        lattice_vector.e12 = R12;
        lattice_vector.e13 = R13;
        lattice_vector.e21 = R21;
        lattice_vector.e22 = R22;
        lattice_vector.e23 = R23;
        lattice_vector.e31 = R31;
        lattice_vector.e32 = R32;
        lattice_vector.e33 = R33;

        inv_lat = lattice_vector.Inverse();

        ModuleBase::Vector3<double> direct_vec, cartesian_vec;
        cartesian_vec.x = cx;
        cartesian_vec.y = cy;
        cartesian_vec.z = cz;
        direct_vec = cartesian_vec * inv_lat;
        dx = direct_vec.x;
        dy = direct_vec.y;
        dz = direct_vec.z;
        return;
    }

    template<typename T>
    static ModuleBase::Vector3<T> latvec_projection(const std::array<ModuleBase::Vector3<T>,3> &latvec)
    {
        ModuleBase::Vector3<T> proj;
        proj.x = std::abs( latvec[0] * (latvec[1] ^ latvec[2]).normalize() );
        proj.y = std::abs( latvec[1] * (latvec[2] ^ latvec[0]).normalize() );
        proj.z = std::abs( latvec[2] * (latvec[0] ^ latvec[1]).normalize() );
        return proj;
    } 
};

} // namespace ModuleBase

#endif

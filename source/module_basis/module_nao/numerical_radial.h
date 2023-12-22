#ifndef NUMERICAL_RADIAL_H_
#define NUMERICAL_RADIAL_H_

#include <cassert>
#include <string>
#include <memory>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"

/**
 * @brief A class that represents a numerical radial function.
 *
 * This class is designed to be the container for the radial part of
 * numerical atomic orbitals, Kleinman-Bylander beta functions, and all
 * other similar numerical radial functions in three-dimensional space,
 * each of which is associated with some angular momentum l and whose r
 * and k space values are related by an l-th order spherical Bessel transform.
 *
 * A NumericalRadial object can be initialized by "build", which requires
 * the angular momentum, the number of grid points, the grid and the
 * corresponding values. Grid does not have to be uniform. One can initialize
 * the object in either r or k space. After initialization, one can set the
 * grid in the other space via set_grid or set_uniform_grid. Values in the
 * other space are automatically computed by a spherical Bessel transform.
 *
 * Usage:
 *
 *     int l = 1;
 *     int itype = 3;
 *     int izeta = 5;
 *     std::string symbol = "Au";
 *
 *     // Prepares the grid & values to initialize the objects
 *     int sz = 2000;
 *     double dr = 0.01;
 *     double* grid = new double[sz];
 *     for (int ir = 0; ir != sz; ++ir) {
 *         grid[ir] = ir * dr; 
 *         f[ir] = std::exp(-grid[ir] * grid[ir]);
 *     }
 *     // grid does not necessarily have to be uniform; it just
 *     // has to be positive and strictly increasing.
 *
 *     // The class will interpret the input values as r^p * F(r)
 *     // where F is the underlying radial function that the class object
 *     // actually represents.
 *     int p1 = 0;
 *     int p2 = -2;
 *
 *     NumericalRadial chi1, chi2;
 *     chi1.build(0, true, sz, grid, f, p1);
 *     chi2.build(2, true, sz, grid, f, p2);
 *
 *     // Now chi1 represents exp(-r^2); chi2 actually represents
 *     // r^2*exp(-r^2), even though the values stored is also exp(-r^2).
 *
 *     // Adds the k-space grid.
 *     chi1.set_uniform_grid(false, sz, PI/dr, 't');
 *     chi2.set_uniform_grid(false, sz, PI/dr, 't');
 *     // k-space values are automatically computed above
 *
 */
class NumericalRadial
{
public:
    NumericalRadial() = default;
    NumericalRadial(NumericalRadial const&); ///< Deep-copy grid & values

    /// Deep-copy grid & values
    NumericalRadial& operator=(NumericalRadial const&);

    ~NumericalRadial();

    /**
     * @brief Initializes the object by providing the grid & values in one space.
     *
     * @param[in] l             Angular momentum.
     * @param[in] for_r_space   Specifies whether the input corresponds to r or k space.
     * @param[in] ngrid         Number of grid points.
     * @param[in] grid          Grid points, must be positive & strictly increasing.
     * @param[in] value         Values on the grid.
     * @param[in] p             Implicit exponent in input values, see @ref pr_ & @ref pk_.
     * @param[in] izeta         Multiplicity index of radial functions of the same itype and l.
     * @param[in] symbol        Chemical symbol.
     * @param[in] itype         Index for the element in calculation.
     * @param[in] init_sbt      If true, internal SphericalBesselTransformer will be initialized.
     *
     * @note init_sbt is only useful when the internal SphericalBesselTransformer (sbt_) is
     *       null-initialized; The function will NOT reset sbt_ if it is already usable.
     */
    void build(const int l,
               const bool for_r_space,
               const int ngrid,
               const double* const grid,
               const double* const value,
               const int p = 0,
               const int izeta = 0,
               const std::string symbol = "",
               const int itype = 0,
               const bool init_sbt = true
    );

    /**
     * @brief Overwrites the content of a Numerical_Orbital_Lm object with the current object.
     *
     * This function provides an interface to the corresponding object in module_ao.
     */
    void to_numerical_orbital_lm(Numerical_Orbital_Lm& orbital_lm, 
                                 const double lcao_ecut = 2000.0,
                                 const double lcao_dk = 0.01) const;

    /** 
     * @brief Sets a SphericalBesselTransformer.
     * 
     * By default the class uses an internal SphericalBesselTransformer, but one can optionally
     * use a shared one. This could be beneficial when there are a lot of NumericalRadial objects
     * whose grids have the same size.
     *
     * @param[in] sbt       An external transformer.
     * @param[in] update    Specifies whether and how values are recomputed with the new transformer.
     *                      Accepted values are:
     *                      *  0: does not recompute values;
     *                      *  1: calls a forward transform;
     *                      * -1: calls a backward transform.
     */
    void set_transformer(ModuleBase::SphericalBesselTransformer sbt, int update = 0);

    /**
     * @brief Sets up a grid.
     *
     * This function can be used to set up the grid which is absent in "build" (in which
     * case values on the new grid are automatically computed by a spherical Bessel transform)
     * or interpolate on an existing grid to a new grid.
     *
     * @param[in] for_r_space   Specifies whether to set grid for the r or k space.
     * @param[in] ngrid         Number of grid points.
     * @param[in] grid          Grid points, must be positive & strictly increasing.
     * @param[in] mode          Specifies how values are updated, could be 'i' or 't':
     *                          - 'i': New values are obtained by interpolating and zero-padding
     *                                 the existing values from current space. With this option,
     *                                 it is an error if the designated space does not have a grid;
     *                          - 't': New values are obtained via transform from the other space.
     *                                 With this option, it is an error if the other space does not
     *                                 have a grid.
     */
    void set_grid(const bool for_r_space, const int ngrid, const double* const grid, const char mode = 'i');

    /**
     * @brief Sets up a uniform grid.
     *
     * The functionality of this function is similar to @ref set_grid, except that
     * the new grid is a uniform grid specified by the cutoff and the number of grid
     * points given by
     *
     *                    cutoff
     *      grid[i] = i * -------
     *                    ngrid-1
     *
     * @see set_grid
     *
     * If enable_fft is true, this function will not only set up the grid & values
     * in the designated space, but also sets the grid in the other space such that
     * the r & k grids are FFT-compliant (and updates values via a FFT-based spherical
     * Bessel transform).
     */
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);

    /**
     * @brief Updates values on an existing grid.
     *
     * This function does not alter the grid; it merely updates values on the existing
     * grid. The number of values to read from "value" is nr_ or nk_. Values of the
     * other space will be automatically updated if they exist.
     *
     * @note This function does not check the index bound; use with care!
     */
    void set_value(const bool for_r_space,
                   const double* const value,
                   const int p
    );

    /// Removes the grid & values in r or k space.
    void wipe(const bool r_space = true, const bool k_space = true);

    //! Computes the radial table for two-center integrals.
    /*!
     * TODO This function shall be removed from this class in the future;
     *     its functionality should be moved to TwoCenterTable class.
     *
     *  This function requires that both "this" and "ket" have existing kgrid_, and the
     *  two kgrid_ be identical.
     *
     *  On finish, table is filled with values on tabgrid. If tabgrid is not given, the
     *  rgrid_ of "this" is assumed (it would be an error if rgrid_ does not exist in
     *  this case).
     *
     *  op could be:
     *
     *  - 'S' or 'I': overlap integral
     *
     *                          / +inf     2
     *              table[ir] = |      dk k  f(k) g(k) j (k*r[ir])
     *                          /  0                    l
     *
     *  - 'T': kinetic integral table
     *
     *                          / +inf     4
     *              table[ir] = |      dk k  f(k) g(k) j (k*r[ir])
     *                          /  0                    l
     *
     *  - 'U': Coulomb integral table. This is slightly different from overlap or
     *         kinetic integral that in this case the two-center integral is a
     *         double integral:
     *
     *                  /        f(r) g(r'-R)
     *                  | dr dr' ------------
     *                  /          |r-r'|
     *
     *         The corresponding table is
     *
     *                          / +inf
     *              table[ir] = |      dk  f(k) g(k) j (k*r[ir])
     *                          /  0                  l
     *
     *                                                                                  */
    void radtab(const char op,              //!< [in] operator, could be:
                                            //!<        - 'S' or 'I': overlap
                                            //!<        - 'T': kinetic
                                            //!<        - 'U': Coulomb
                const NumericalRadial& ket, //!< [in] the other NumericalRadial object with which
                                            //!       the two-center integral is computed
                const int l,                //!< [in] angular momentum of the table
                double* const table,        //!< [out] on finish, contain the computed table
                const int nr_tab,           //!< [in] size of table grid
                const double rmax_tab,      //!< [in] cutoff radius of table grid
                const bool deriv = false    //!< [in] if true, calculates the derivative of the table
    ) const;

    /**
     * @brief Normalizes the radial function.
     *
     * The radial function is normalized such that
     *
     *      / +inf     2
     *      |      dx x  f(x) = 1
     *      /  0
     *
     * where x is r or k. The integral is evaluated with Simpson's rule. Values in the other space
     * are updated automatically via a spherical Bessel transform.
     */
    void normalize(bool for_r_space = true);

    /**
     * @name Getters
     */
    ///@{
    std::string const& symbol() const { return symbol_; }
    int itype() const { return itype_; }
    int izeta() const { return izeta_; }
    int l() const { return l_; }
    int nr() const { return nr_; } // paired with rmax(), not rcut!
    int nk() const { return nk_; }
    double rcut() const { return rgrid_[std::min(ircut_, nr_-1)]; } ///< padded zeros ignored
    double kcut() const { return kgrid_[std::min(ikcut_, nk_-1)]; }
    double rmax() const { return rgrid_[nr_-1]; } ///< padded zeros considered
    double kmax() const { return kgrid_[nk_-1]; }
    const double* rgrid() const { return rgrid_; }
    const double* kgrid() const { return kgrid_; }
    const double* rvalue() const { return rvalue_; }
    const double* kvalue() const { return kvalue_; }
    double pr() const { return pr_; }
    double pk() const { return pk_; }
    bool is_fft_compliant() const { return is_fft_compliant_; }
    ModuleBase::SphericalBesselTransformer sbt() const { return sbt_; }

    double rgrid(int ir) const { return rgrid_[ir]; }
    double kgrid(int ik) const { return kgrid_[ik]; }
    double rvalue(int ir) const { return rvalue_[ir]; }
    double kvalue(int ik) const { return kvalue_[ik]; }
    ///@}

private:
    std::string symbol_ = "";   ///< chemical symbol
    int itype_ = 0;             ///< element index in calculation
    int l_ = -1;                ///< angular momentum
    int izeta_ = 0;             ///< further index for NumericalRadial objects with the same itype_and l_

    int nr_ = 0;                ///< number of r-space grid points
    int nk_ = 0;                ///< number of k-space grid points

    double* rgrid_ = nullptr;   ///< r-space grid
    double* kgrid_ = nullptr;   ///< k-space grid

    /**
     * @brief Index of the first trailing zero.
     *
     * A numerical radial function might be zero-padded for the sake of
     * FFT-based spherical Bessel transform algorithms. The following two
     * variables keep track of the actual cutoff radius. Specifically,
     * if there are no trailing zeros in rvalues_, then ircut_ = nr_;
     * if there are trailing zeros, then ircut_ is the index of the first
     * trailing zero. For example, 
     * rvalues_ = {1, 2, 3, 0, 0, 0} -> ircut_ = 3
     * rvalues_ = {1, 2, 3, 4, 5, 6} -> ircut_ = 6
     * rvalues_ = {0, 0, 0, 0, 0, 0} -> ircut_ = 0
     */
    int ircut_ = 0;
    int ikcut_ = 0;

    double* rvalue_ = nullptr;  ///< r-space value
    double* kvalue_ = nullptr;  ///< k-space value

    /**
     * @brief A flag that tells whether the r & k grids are FFT-compliant.
     *
     * r & k grids are considered FFT-compliant if they
     * -# have the same number of grid points;
     * -# are both uniform;
     * -# both starts from 0;
     * -# satisfy dr*dk = pi/(N-1) where N >= 2 is the number of each grid points.
     *
     * If the grids are FFT-compliant, spherical Bessel transforms are performed
     * with an FFT-based algorithm. Otherwise, the transforms are performed with
     * numerical integration (Simpson's rule).
     */
    bool is_fft_compliant_ = false;

    /**
     * @name Implicit exponents in values
     *
     * Sometimes a radial function is given in the form of pow(r,p) * F(r) rather
     * than F(r) (same applies to k). For example, the Kleinman-Bylander beta
     * functions are often given as r*beta(r) instead of bare beta(r). Very often
     * using r*beta(r) is adequate; there's no need to get bare beta(r) at all.
     *
     * This class takes care of this situation. When building the object, one can
     * optionally provide an exponent p so that the input values are interpreted as
     * pow(x[i],p) * F(x[i]), where F(x) is what the object actually represents.
     * pr_ & pk_ keep track of these exponents within r & k values. They are taken
     * taken account automatically during spherical Bessel transforms.
     */
    ///@{
    /**
     * This parameter affects how this class interprets rvalues_. Specifically,
     * rvalues_[ir] is interpreted as pow(rgrid_[ir], pr_) * F(rgrid_[ir])
     * during spherical Bessel transforms.
     */
    int pr_ = 0; ///< implicit exponent in rvalues_

    /**
     * This parameter affects how this class interprets kvalues_. Specifically,
     * kvalues_[ik] is interpreted as pow(kgrid_[ik], pk_) * F(kgrid_[ik])
     * during spherical Bessel transforms.
     */
    int pk_ = 0; ///< implicit exponent in kvalues_
    ///@}

    /// An object that provides spherical Bessel transforms
    ModuleBase::SphericalBesselTransformer sbt_{nullptr};

    /**
     * @brief Transforms the r-space values to get k-space values, or vice versa.
     *
     * The grid & values where the transform is initiated must exist; this function
     * does nothing if grid in the destination space does not exist.
     *
     * forward : r to k
     * backward: k to r
     */
    void transform(const bool forward);

    /// Updates ircut_ or ikcut_.
    void set_icut(const bool for_r_space, const bool for_k_space, const double tol = 1e-15);

    /// Checks whether a grid is uniform.
    bool is_uniform(const int n, const double* const grid, const double tol) const;

    /**
     * @brief Checks whether the given two grids are FFT-compliant.
     *
     * @see is_fft_compliant_
     */
    bool is_fft_compliant(const int nr, const double* const rgrid, const int nk, const double* const kgrid) const;
};

#endif

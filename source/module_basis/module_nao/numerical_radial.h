#ifndef NUMERICAL_RADIAL_H_
#define NUMERICAL_RADIAL_H_

#include <cassert>
#include <string>

#include "module_base/spherical_bessel_transformer.h"

//! A class that represents a numerical radial function.
/*!
 *  This class is supposed to be the underlying container for numerical
 *  atomic orbitals, Kleinman-Bylander beta functions, and all other
 *  similar numerical radial functions in three-dimensional space which
 *  is associated with some angular momentum l and whose r & k space
 *  values are related by an l-th order spherical Bessel transform.
 *
 *  A NumericalRadial object can be initialized by "build", which requires
 *  the angular momentum, the number of grid points, the grid and the
 *  corresponding values. One can initialize the object in either r or
 *  k space. After initialization, one can set the grid in the other
 *  space. Values in the other space are automatically computed by a
 *  spherical Bessel transform.
 *
 *  This class also provides a convenient interface to compute the radial
 *  table for various two-center intergrals.
 *
 *  Usage:
 *
 *      int l = 1;
 *      int itype = 3;
 *      int izeta = 5;
 *      std::string symbol = "Au";
 *
 *      // Prepares the grid & values to initialize the objects
 *      int sz = 2000;
 *      double dr = 0.01;
 *      double* grid = new double[sz];
 *      for (int ir = 0; ir != sz; ++ir) {
 *          double r = ir * dr;
 *          grid[ir] = r;
 *          f[ir] = std::exp(-r*r);
 *      }
 *
 *      // The class will interpret the input values as r^p * F(r)
 *      // where F is the underlying radial function that the class object
 *      // actually represents.
 *      int p1 = 0;
 *      int p2 = -2;
 *
 *      NumericalRadial chi1, chi2;
 *      chi1.build(0, true, sz, grid, f, p1);
 *      chi2.build(2, true, sz, grid, f, p2);
 *
 *      // Now chi1 represents exp(-r^2); chi2 actually represents
 *      // r^2*exp(-r^2), even though the values stored is also exp(-r^2).
 *
 *      // Adds the k-space grid.
 *      chi1.set_uniform_grid(false, sz, PI/dr, 't');
 *      chi2.set_uniform_grid(false, sz, PI/dr, 't');
 *      // k-space values are automatically computed above
 *
 *      // calculates various radial tables between chi1 & chi2
 *      double* table = new double[sz];
 *      chi1.radtab('S', chi2, 0, table);
 *
 *                                                                          */
class NumericalRadial
{

  public:
    NumericalRadial();
    NumericalRadial(NumericalRadial const&); //!< deep copy

    //! deep copy
    NumericalRadial& operator=(NumericalRadial const&);

    ~NumericalRadial();

    //! Initializes the object by providing the grid & values in one space.
    void build(const int l,               //!< angular momentum
               const bool for_r_space,    //!< specifies whether the input corresponds to r or k space
               const int ngrid,           //!< number of input grid points
               const double* const grid,  //!< must be strictly increasing, and every element must be larger than zero
               const double* const value, //!< values on the grid
               const int p = 0,           //!< exponent of the implicit power term in input values, @see @ref group1
               const int izeta = 0,       //!< index for the multiplicity of radial functions of the same itype and l
               const std::string symbol = "", //!< chemical symbol
               const int itype = 0            //!< index for the element in calculation
    );

    //! Sets a SphericalBesselTransformer.
    /*!
     *  By default the class uses an internal SphericalBesselTransformer,
     *  but one can optionally use an external one. This could be beneficial
     *  when there are a lot of NumericalRadial objects whose grids are all
     *  FFT-compliant and have the same size. In that case, one can set up
     *  an external transformer, set the FFTW plan flag to FFTW_MEASURE,
     *  and have all NumericalRadial objects use the external transformer.
     *                                                                      */
    void set_transformer(ModuleBase::SphericalBesselTransformer* sbt = nullptr, //!< pointer to external transformer.
                                                                                //!< nullptr instructs the object to use
                                                                                //!< an an internal transformer.
                         int update = 0 //!< Specifies whether and how values are recomputed with the new transformer.
                                        //!< Accepted values are:
                                        //!< *  0: does not recompute values;
                                        //!< *  1: calls a forward transform
                                        //!< * -1: calls a backward transform
    );

    //! Sets up a new grid
    void set_grid(const bool for_r_space,   //!< specifies whether to set grid for the r or k space
                  const int ngrid,          //!< number of grid points
                  const double* const grid, //!< must be stricly increasing, and every element must be larger than zero
                  const char mode = 'i'     //!< specifies how values are updated, could be 'i' or 't'.
                                            //!< - 'i': new values are obtained by interpolating and zero-padding
                                            //!<        the existing values from current space.
                                            //!< - 't': new values are obtained via transform from the other space
    );

    /*!
     *  Sets up a new uniform grid by
     *
     *                    cutoff
     *      grid[i] = i * -------
     *                    ngrid-1
     *
     *  @see set_grid
     *
     *  If enable_fft is true, this function will not only set up the grid & values
     *  in the designated space, but also sets the grid (and updates values accordingly)
     *  in the other space such that r & k grids are FFT-compliant.
     *                                                                                  */
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);

    //! Updates values on an existing grid.
    /*!
     *  The number of values to read from "value" is nr_ or nk_ depending on for_r_space.
     *  Values of the other space will also be updated if they exist.
     *  This function does not check the index bound; use with care!
     *                                                                                  */
    void set_value(const bool for_r_space,    //!< specifies whether to set grid for the r or k space
                   const double* const value, //!< new values
                   const int p                //!< see @ref group1
    );

    //! Removes the grid & values from one space.
    void wipe(const bool r_space /*! specifies whether to wipe off the r or k space info */);

    //! Saves the data to file (what data, in what format?)
    void save(std::string file = "" /*! file name */) const;

    //! Computes the radial table for two-center integrals.
    /*!
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
    void radtab(const char op,                 //!< [in] operator, could be:
                                               //!<        - 'S' or 'I': overlap
                                               //!<        - 'T': kinetic
                                               //!<        - 'U': Coulomb
                const NumericalRadial& ket,    //!< [in] the other NumericalRadial object with which
                                               //!       the two-center integral is computed
                const int l,                   //!< [in] angular momentum of the table
                double* const table,           //!< [out] on finish, contain the computed table
                const int nr_tab,              //!< [in] size of table grid
                const double* const rgrid_tab, //!< [in] grid on which the table is calculated.
                const bool deriv = false       //!< [in] if true, calculates the derivative of the table
    ) const;

    //! Normalizes the radial function.
    /*!
     *  The radial function is normalized such that
     *
     *      / +inf     2
     *      |      dx x  f(x) = 1
     *      /  0
     *
     *  where x is r or k.
     *                                                                                  */
    void normalize(bool for_r_space = true);

    /*!
     *  @name Getters
     *                                                                                  */
    ///@{
    //! gets symbol_
    std::string const& symbol() const { return symbol_; }

    //! gets itype_
    int itype() const { return itype_; }

    //! gets izeta_
    int izeta() const { return izeta_; }

    //! gets the angular momentum
    int l() const { return l_; }

    //! gets the number of r-space grid points
    int nr() const { return nr_; }

    //! gets the number of k-space grid points
    int nk() const { return nk_; }

    //! gets r-space grid cutoff distance
    double rcut() const { return rgrid_ ? rgrid_[nr_ - 1] : 0.0; }

    //! gets k-space grid cutoff distance
    double kcut() const { return kgrid_ ? kgrid_[nk_ - 1] : 0.0; }

    //! gets the pointer to r-space grid points
    const double* ptr_rgrid() const { return rgrid_; }

    //! gets the pointer to k-space grid points
    const double* ptr_kgrid() const { return kgrid_; }

    //! gets the pointer to r-space values
    const double* ptr_rvalue() const { return rvalue_; }

    //! gets the pointer to k-space values
    const double* ptr_kvalue() const { return kvalue_; }

    //! gets the exponent of the pre-multiplied power term in rvalues_. @see pr_
    double pr() const { return pr_; }

    //! gets the exponent of the pre-multiplied power term in kvalues_. @see pk_
    double pk() const { return pk_; }

    //! gets the flag for FFT-compliancy. @see is_fft_compliant_
    bool is_fft_compliant() const { return is_fft_compliant_; }

    //! gets the pointer to the SphericalBesselTransformer
    const ModuleBase::SphericalBesselTransformer* ptr_sbt() const { return sbt_; }

    ///@}

  private:
    std::string symbol_ = ""; //!< chemical symbol
    int itype_ = 0;           //!< element index in calculation
    int l_ = -1;              //!< angular momentum
    int izeta_ = 0;           //!< further index for NumericalRadial objects with the same itype_and l_

    int nr_ = 0; //!< number of r-space grid points
    int nk_ = 0; //!< number of k-space grid points

    double* rgrid_ = nullptr; //!< r-space grid
    double* kgrid_ = nullptr; //!< k-space grid

    //! A flag that tells whether the r & k grids are FFT-compliant.
    /*!
     *  r & k grids are considered FFT-compliant if they
     *  1. have the same number of grid points;
     *  2. are both uniform;
     *  3. both starts from 0;
     *  4. satisfy dr*dk = pi/(N-1) where N >= 2 is the number of each grid points
     *                                                                              */
    bool is_fft_compliant_ = false;

    double* rvalue_ = nullptr; //!< r-space value
    double* kvalue_ = nullptr; //!< k-space value

    /*!
     *  @name Exponents of the implicit power terms
     *
     *  Sometimes a radial function is given in the form of pow(r,p) * F(r) rather
     *  than F(r) (same applies to k). For example, the Kleinman-Bylander beta
     *  functions are often given as r*beta(r) instead of bare beta(r). Very often
     *  r*beta(r) is adequate; bare beta(r) is not necessary at all.
     *
     *  This class takes care of this situation. When building the object, one can
     *  optionally provide an exponent p so that the values are interpreted as
     *  pow(r[i],p) * F(r[i]). pr_ & pk_ keep track of these exponents within r & k
     *  values. They are automatically taken account during spherical Bessel
     *  transforms.
     *                                                                              */
    ///@{
    /*!
     *  This parameter affects how this class interprets rvalues_. Specifically,
     *  rvalues_[ir] is interpreted as pow(rgrid_[ir], pr_) * F(rgrid_[ir])
     *  during spherical Bessel transforms.
     *                                                                              */
    int pr_ = 0; //!< exponent of the implicit power term in rvalues_

    /*!
     *  This parameter affects how this class interprets kvalues_. Specifically,
     *  kvalues_[ik] is interpreted as pow(kgrid_[ik], pk_) * F(kgrid_[ik])
     *  during spherical Bessel transforms.
     *                                                                              */
    int pk_ = 0; //!< exponent of the implicit power term in kvalues_
    ///@}

    //! Pointer to the object that provides spherical Bessel transforms
    /*!
     *  @see set_transformer
     *                                                                              */
    ModuleBase::SphericalBesselTransformer* sbt_;

    //! A flag that marks the ownership of sbt_
    bool use_internal_transformer_;

    //! Transforms the r-space values to get k-space values, or vice versa.
    /*!
     *  The grid & values where the transform is initiated must exist; this function
     *  does nothing if grid in the destination space does not exist.
     *
     *  forward : r to k
     *  backward: k to r
     *                                                                              */
    void transform(const bool forward);

    //! Checks whether the given two grids are FFT-compliant
    /*!
     *  @see is_fft_compliant
     *                                                                              */
    bool is_fft_compliant(const int nr, const double* const rgrid, const int nk, const double* const kgrid) const;
};

#endif

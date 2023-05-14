#include <string>
#include <cassert>
#include "module_base/spherical_bessel_transformer.h"

class NumericalRadial {

public:
   
    NumericalRadial();
    NumericalRadial(NumericalRadial const&);

    //! Performs a deep copy
    NumericalRadial& operator=(NumericalRadial const&);

    ~NumericalRadial();

    //! Initializes the object by providing the grid & values in one space.
    void build(
            const int l,                    //!< [in] angular momentum
            const char r_or_k,              //!< [in] specifies whether the input corresponds to r or k space
            const int ngrid,                //!< [in] number of input grid / values points
            const double* const grid,       //!< [in] grid
            const double* const value,      //!< [in] values on the grid
            const int p = 0,                //!< [in] exponent of the implicit power term in input values, @see @ref group1 
            const int itype = 0,            //!< [in] usually the index for elements
            const int ichi = 0,             //!< [in] further index after itype and l
            const std::string symbol = ""   //!< [in] usually the chemical symbol
    );

    //! Initializes the object by providing the grids in both space and values in one space.
    void build(
            const int l,                    //!< [in] angular momentum
            const int nr,                   //!< [in] number of r-space grid points
            const double* const rgrid,      //!< [in] r-space grid
            const int nk,                   //!< [in] number of k-space grid points
            const double* const kgrid,      //!< [in] k-space grid
            const char r_or_k,              //!< [in] specifies whether the values corresponds to r or k space
            const double* const value,      //!< [in] values on the grid
            const int p = 0,                //!< [in] exponent of the implicit power term in input values, @see @ref group1
            const int itype = 0,            //!< [in] usually the index for elements
            const int ichi = 0,             //!< [in] further index after itype and l
            const std::string symbol = ""   //!< [in] usually the chemical symbol
    );

    //! Sets a SphericalBesselTransformer.
    /*!
     *  By default the class uses an internal SphericalBesselTransformer,
     *  but one can optionally use an external one. This could be beneficial
     *  when there are a lot of NumericalRadial objects whose grids are all
     *  FFT-compliant and have the same size. In that case, one can set up 
     *  an external transformer, set the FFTW plan flag to FFTW_MEASURE, 
     *  and have all NumericalRadial objects use this transformer.
     *
     *  If sbt is nullptr, the class will use an internal one.
     *                                                                      */
    void set_transformer(
            ModuleBase::SphericalBesselTransformer* sbt = nullptr, //!< pointer to external transformer
            int update = 0 //!< specifies whether and how values are recomputed with the new transformer
    );

    //! Sets up a new grid
    void set_grid(
            const char r_or_k,          //!< [in] 'r' or 'k'
            const int ngrid,            //!< [in] number of grid points
            const double* const grid,   //!< [in] grid
            const char mode = 'i'       //!< [in] 'i' or 't'.
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
     *  If enable_fft is true, this function will first set up the grid & values
     *  in r_or_k space, and then sets the FFT-compliant grid in the other space,
     *  and transform to get new values.
     *                                                                                  */
    void set_grid(
            const char r_or_k, 
            const int ngrid, 
            const double cutoff,
            const char mode = 'i',
            const bool enable_fft = false
    );

    //! Updates values on an existing grid.
    /*!
     *  Values of the other space will also be updated if it exist.
     *                                                                                  */
    void set_value(
            const char r_or_k,
            const double* const value
    );

    //! Removes the grid & values from one space.
    void wipe(const char r_or_k);

    //! Saves the data to file (what data, in what format?)
    void save(const std::string & file = "", const bool include_header = true) const;

    //! Computes the radial table for the two-center integral.
    /*!
     *  Currently this function requires that "this" and "ket" have exactly the same 
     *  grid and are FFT-compliant. On finish, table is filled with values on the same 
     *  rgrid_ of each object.
     *
     *  op could be:
     *
     *  - 'S' or 'I': overlap integral
     *
     *          
     *          /
     *          | f(r) g(r-R) dr
     *          /
     *
     *                      / +inf     2
     *          table[ir] = |      dk k  f(k) g(k) j (k*r[ir])
     *                      /  0                    l
     *
     *  - 'T': kinetic integral table
     *
     *                      / +inf     4
     *          table[ir] = |      dk k  f(k) g(k) j (k*r[ir])
     *                      /  0                    l
     *
     *  - 'U': Coulomb integral table. This is slightly different from overlap or 
     *         kinetic integral that in this case the two-center integral is a 
     *         double integral:
     *
     *         /        f(r) g(r'-R)
     *         | dr dr' ------------
     *         /          |r-r'|
     *
     *         The corresponding table is
     *
     *                      / +inf    
     *          table[ir] = |      dk  f(k) g(k) j (k*r[ir])
     *                      /  0                    l
     *
     *                                                                                  */
    void radtab(
            const char op,
            const NumericalRadial& ket,
            const int l, 
            double* const table,
            const bool deriv = false
    );

    /*! 
     *  @name Getters
     *                                                                                  */
    ///@{
    //! gets symbol_
    std::string const& symbol() const { return symbol_; }

    //! gets the angular momentum
    int l() const { return l_; }

    //! gets the number of r-space grid points
    int nr() const { return nr_; }

    //! gets the number of k-space grid points
    int nk() const { return nk_; }

    //! gets r-space grid cutoff distance
    double rcut() const { return rgrid_[nr_-1]; }

    //! gets k-space grid cutoff distance
    double kcut() const { return kgrid_[nk_-1]; }

    //! gets the pointer to r-space grid points
    double* ptr_rgrid() const { return rgrid_; }

    //! gets the pointer to k-space grid points
    double* ptr_kgrid() const { return kgrid_; }

    //! gets the pointer to r-space values
    double* ptr_rvalue() const { return rvalue_; }

    //! gets the pointer to k-space values
    double* ptr_kvalue() const { return kvalue_; }

    //! gets the ir-th r-space grid point
    double rgrid(const int ir) const { assert(ir<nr_); return rgrid_[ir]; }

    //! gets the ik-th k-space grid point
    double kgrid(const int ik) const { assert(ik<nk_); return kgrid_[ik]; }

    //! gets the value on the ir-th r-space grid point
    double rvalue(const int ir) const { assert(ir<nr_); return rvalue_[ir]; }

    //! gets the value on the ik-th k-space grid point
    double kvalue(const int ik) const { assert(ik<nk_); return kvalue_[ik]; }

    //! gets the exponent of the pre-multiplied power term in rvalues_. @see pr_
    double pr() const { return pr_; }

    //! gets the exponent of the pre-multiplied power term in kvalues_. @see pk_
    double pk() const { return pk_; }
    ///@}

private:

    int l_ = -1; //!< angular momentum

    int nr_ = 0; //!< number of r-space grid points
    int nk_ = 0; //!< number of k-space grid points

    double* rgrid_ = nullptr; //!< r-space grid
    double* kgrid_ = nullptr; //!< k-space grid

    double* rvalue_ = nullptr; //!< r-space value
    double* kvalue_ = nullptr; //!< k-space value

    //! A flag that tells whether the r & k grids are FFT-compliant.
    /*!
     *  r & k grids are considered FFT-compliant if they
     *  1. have the same number of grid points;
     *  2. are both uniform;
     *  3. both starts from 0;
     *  4. satisfy dr*dk = pi/(N-1) where N >= 2 is the number of each grid points 
     *                                                                              */
    bool is_fft_compliant_ = false;

    //! An object that provides spherical Bessel transform
    /*! 
     *  The SphericalBesselTransformer class is designed to efficiently perform 
     *  many transforms of the same size.
     *                                                                              */
    ModuleBase::SphericalBesselTransformer* sbt_ = nullptr;

    //! A flag that tells the ownership of sbt_
    bool use_internal_transformer_ = false;

    /*! 
     *  @defgroup group1 Exponents of the implicit power terms
     *
     *  Sometimes a radial function is given in the form of pow(r,p) * F(r) rather
     *  than F(r) (same applies to k). For example, the Kleinman-Bylander beta 
     *  functions are often given as r*beta(r) instead of bare beta(r), and very 
     *  often all one needs is r*beta(r) & beta(k); one never needs the bare beta(r).
     *
     *  This class takes care of this situation. When building the object, one can 
     *  specify the exponent p and just pass pow(r[i],p) * F(r[i]) (or the k-space 
     *  counterpart) to the value. pr_ & pk_ keep track of these exponents within
     *  r & k values. They are automatically taken account during spherical Bessel 
     *  transforms. 
     *                                                                              */
    ///@{
    /*! Interprets rvalues_[ir] as pow(rgrid_[ir], pr_) * F(rgrid_[ir]) */
    int pr_ = 0; //!< exponent of the implicit power term in rvalues_

    /*! Interprets kvalues_[ik] as pow(kgrid_[ik], pk_) * F(kgrid_[ik]) */
    int pk_ = 0; //!< exponent of the implicit power term in kvalues_
    ///@}

    std::string symbol_ = ""; //!< usually the element symbol

    int itype_ = 0; //!< usually the index for element

    int ichi_ = 0; //!< further index for NumericalRadial objects with the same itype_ and l_

    //! Applies a spherical Bessel transform to r(k) space values to get k(r) space values.
    /*! 
     *  r & k grids must exist; output value array must be pre-allocated.
     *
     *  forward : r to k
     *  backward: k to r
     *                                                                              */
    void transform(const bool forward);

    //! Checks whether r & k grids are FFT-compliant and set the corresponding flag.
    /*!
     *  @see is_fft_compliant
     *                                                                              */
    void check_fft_compliancy();
};



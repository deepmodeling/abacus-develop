#ifndef RADIAL_SET_H_
#define RADIAL_SET_H_

#include <map>
#include <string>
#include <utility>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/numerical_radial.h"
#include "module_basis/module_ao/ORB_nonlocal.h"

//! An abstract class representing a related set of numerical radial functions.
/*!
 *  This abstract class represents a set of numerical radial functions from
 *  a single file. This is supposed to be the base class for concrete classes
 *  like AtomicRadials and BetaRadials, in which case they represent the set
 *  of all radial functions of the numerical atomic orbitals and Kleinman-
 *  Bylander beta functions of a single element respectively.
 *
 *  @see AtomicRadials BetaRadials
 *                                                                          */
class RadialSet
{
  public:
    RadialSet();
    RadialSet(const RadialSet&);            //!< deep copy
    RadialSet& operator=(const RadialSet&); //!< deep copy
    virtual RadialSet* clone() const = 0;   //!< for polymorphic copy

    virtual ~RadialSet();

    /*! @brief Build the object from a file.
     *
     * Currently only AtomicRadials objects are supposed to used this
     * interface.
     *                                                                      */
    virtual void build(const std::string&,             ///< file name
                       const int = 0,                  ///< the element index in calculation
                       std::ofstream* const = nullptr, ///< output file stream for logging
                       const int = 0                   ///< MPI rank
                       ) {}

    /*! @brief Build from a Numerical_Nonlocal object.
     *
     * This function is supposed to be used by BetaRadials ONLY.
     *                                                                      */
    virtual void build(const Numerical_Nonlocal&,     ///< Numerical_Nonlocal object
                       const int = 0,                 ///< the element index in calculation
                       std::ofstream* const = nullptr ///< output file stream for logging
                       ) {}

    /*! @name Getters
     *
     *  Get access to private members.
     *                                                                      */
    //!@{
    const std::string& symbol() const { return symbol_; }
    int itype() const { return itype_; }
    int lmax() const { return lmax_; }
    double rcut_max() const;

    int nzeta(const int l) const { return (l >= 0 && l <= lmax_) ? nzeta_[l] : 0; }
    int nzeta_max() const { return nzeta_max_; }
    int nchi() const { return nchi_; }

    const NumericalRadial& chi(const int l, const int izeta);
    const NumericalRadial* cbegin() const { return chi_; }
    const NumericalRadial* cend() const { return chi_ + nchi_; }
    //!@}

    /*! @name property setters for all NumericalRadial objects
     *
     *  @see NumericalRadial
     *                                                                     */
    //!@{
    //! Set a spherical Bessel transformers for all NumericalRadial objects
    //! @see NumericalRadial::set_transformer
    void set_transformer(ModuleBase::SphericalBesselTransformer* const sbt = nullptr, const int update = 0);

    //! Set a common grid for all NumericalRadial objects
    //! @see NumericalRadial::set_grid
    void set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode = 'i');

    //! Set a common uniform grid for all NumericalRadial objects
    //! @see NumericalRadial::set_uniform_grid
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);
    //!@}

  protected:
    std::string symbol_ = ""; //!< usually the chemical symbol
    int itype_ = 0;           //!< usually the index for element in calculation
    int lmax_ = -1;           //!< maximum angular momentum among all NumericalRadial objects

    int* nzeta_ = nullptr; //!< number of NumericalRadial objects for each angular momentum
    int nzeta_max_ = 0;    //!< maximum number of NumericalRadial objects among each angular momentum
    int nchi_ = 0;         //!< total number of NumericalRadial objects

    NumericalRadial* chi_ = nullptr; //!< array of NumericalRadial objects

    int* index_map_ = nullptr; //!< map (l,izeta) to an index in chi_ array

    //! Pointer to the object that provides spherical Bessel transforms
    /*!
     *  All NumericalRadial objects within this class should share the same
     *  spherical Bessel transformer.
     *                                                                      */
    ModuleBase::SphericalBesselTransformer* sbt_;

    //! A flag that marks the ownership of sbt_
    bool use_internal_transformer_;

    //! deallocates memory and reset all class members to default values
    void cleanup();

    //! get the index in chi_ array from (l,izeta)
    int index(const int l, const int izeta) const;

    //! calculate nzeta_max_ and build index_map_ from nzeta_ and lmax_
    void indexing();
};

#endif

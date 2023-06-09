#ifndef RADIAL_SET_H_
#define RADIAL_SET_H_

#include <map>
#include <string>
#include <utility>

#include "module_basis/module_nao/numerical_radial.h"

//! An abstract class representing a related set of numerical radial functions
/*!
 *  This abstract class represents a set of numerical radial functions from
 *  a single file and of the same "itype" number. This is supposed to be the
 *  base class for concrete classes like AtomicRadials and BetaRadials, in
 *  which case "itype" is the element index in calculation, and they represent
 *  the set of all radial functions of the numerical atomic orbitals and
 *  Kleinman-Bylander beta functions of a single element respectively.
 *
 *  @see AtomicRadials BetaRadials
 *                                                                          */
class RadialSet
{
  public:
    RadialSet(){};
    ~RadialSet();

    //! build the set of numerical radial functions from a file
    virtual void build(const std::string& file,                //!< orbital or pseudopotential file
                       const int itype = 0,                    //!< the element index in calculation
                       std::ofstream* const ptr_log = nullptr, //!< output file stream for logging
                       const int rank = 0                      //!< MPI rank
                       )
        = 0;

    /*! @name Getters
     *
     *  Get access to private members.
     *                                                                      */
    //!@{
    const std::string& symbol() const
    {
        return symbol_;
    }
    int itype() const
    {
        return itype_;
    }
    int lmax() const
    {
        return lmax_;
    }

    int nzeta(int l) const
    {
        return nzeta_[l];
    }
    int nzeta_max() const
    {
        return nzeta_max_;
    }
    int nchi() const
    {
        return nchi_;
    }

    double rcut_max() const
    {
        return rcut_max_;
    }

    const NumericalRadial& chi(int l, int izeta);
    //!@}

    /*! @name property setters for all NumericalRadial objects
     *
     *  @see NumericalRadial
     *                                                                     */
    //!@{
    //! Set a spherical Bessel transformers for all NumericalRadial objects
    //! @see NumericalRadial::set_transformer
    void set_transformer(ModuleBase::SphericalBesselTransformer* sbt = nullptr, int update = 0);

    //! Set a common grid for all NumericalRadial objects
    //! @see NumericalRadial::set_grid
    void set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode = 'i');

    //! Set a common uniform grid for all NumericalRadial objects
    //! @see NumericalRadial::set_uniform_grid
    void set_uniform_grid(const bool for_r_space, const int ngrid, const double cutoff, const char mode = 'i');
    //!@}

  protected:
    std::string symbol_ = ""; //!< usually the chemical symbol
    int itype_ = 0;           //!< usually the index for element in calculation
    int lmax_ = -1;           //!< maximum angular momentum among all NumericalRadial objects

    int* nzeta_ = nullptr; //!< number of NumericalRadial objects for each angular momentum
    int nzeta_max_ = 0;    //!< maximum number of NumericalRadial objects among each angular momentum
    int nchi_ = 0;         //!< total number of NumericalRadial objects

    double rcut_max_ = 0; //!< maximum r-space cutoff radius among all NumericalRadial objects

    NumericalRadial* chi_ = nullptr; //!< array of NumericalRadial objects

    int* index_map_ = nullptr; //!< map (l,izeta) to an index in chi_ array

    //! deallocates memory and reset all class members to default values
    void cleanup();

    //! get the index in chi_ array from (l,izeta)
    int index(const int l, const int izeta) const;
};

#endif

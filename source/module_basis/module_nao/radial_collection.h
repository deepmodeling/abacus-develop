#ifndef RADIAL_COLLECTION_H_
#define RADIAL_COLLECTION_H_

#include <numeric>
#include <string>

#include "module_basis/module_nao/radial_set.h"

//! A class that holds the collection of all numerical radial functions of a kind.
/*!
 *  This class is supposed to hold the radial functions of all numerical atomic
 *  orbitals or Kleinman-Bylander beta functions from all elements involved in
 *  calculation. The way it stores the radial functions is type(element)-wise,
 *  i.e. all radial functions of the same itype (element) are stored in a
 *  RadialSet object, and the collection holds an array of such objects.
 *                                                                              */
class RadialCollection
{
  public:
    RadialCollection();
    RadialCollection(const RadialCollection& other); //!< deep copy

    RadialCollection& operator=(const RadialCollection& rhs); //!< deep copy

    ~RadialCollection();

    /// build the collection from (orbital) files
    void build(const int nfile, const std::string* const file, const char type = 'o');

    /// build the collection from Numerical_Nonlocal objects
    void build(const int ntype, Numerical_Nonlocal* const nls);

    /*! @name Getters
     *
     *  Get access to private members (and their properties).
     *                                                                      */
    //!@{
    //! element symbol of a given type
    const std::string& symbol(const int itype) const { return radset_[itype]->symbol(); }

    //! number of RadialSet objects in the collection
    int ntype() const { return ntype_; }

    //! maximum angular momentum of the itype-th RadialSet in the collection
    int lmax(const int itype) const { return radset_[itype]->lmax(); }

    //! maximum angular momentum of all NumericalRadial objects in the collection
    int lmax() const { return lmax_; }

    //! maximum cutoff radius of a give type
    double rcut_max(const int itype) const { return radset_[itype]->rcut_max(); }

    //! maximum cutoff radius of all NumericalRadial objects in the collection
    double rcut_max() const;

    //! number of distinct radial functions of a given type and angular momentum
    int nzeta(const int itype, const int l) const { return radset_[itype]->nzeta(l); }

    //! maximum number of distinct radial functions of a given type among all angular momentum
    int nzeta_max(const int itype) const { return radset_[itype]->nzeta_max(); }

    //! maximum number of distinct radial functions of a given type among all angular momentum
    int nzeta_max() const { return nzeta_max_; }

    //! total number of NumericalRadial objects in the collection
    int nchi() const { return nchi_; }

    //! number of NumericalRadial objects of a given type
    int nchi(const int itype) const { return radset_[itype]->nchi(); }

    //! get access to the NumericalRadial object with given type, angular momentum and zeta number
    const NumericalRadial& operator()(const int itype, const int l, const int izeta) const
    {
        assert(itype >= 0 && itype < ntype_);
        return radset_[itype]->chi(l, izeta);
    }

    //! get access to the RadialSet object with given type
    const RadialSet& operator()(const int itype) const
    {
        assert(itype >= 0 && itype < ntype_);
        return *radset_[itype];
    }
    //!@}

    /*! @name Iterators
     *
     *  Enable iteration through all NumericalRadial objects in the collection.
     *  Objects are sorted by l first, by itype next, by izeta last.
     *                                                                      */
    //!@{
    const NumericalRadial** cbegin() const
    {
        assert(ntype_ > 0);
        return iter_;
    }

    const NumericalRadial** cend() const
    {
        assert(ntype_ > 0);
        return iter_ + nchi_;
    }

    //! *(this->cbegin(l)) returns the address of the first NumericalRadial object with angular momentum l
    const NumericalRadial** cbegin(const int l) const
    {
        assert(ntype_ > 0 && l >= 0 && l <= lmax_);
        return iter_ + std::accumulate(nl_, nl_ + l, 0);
    }

    //! *(this->cbegin(l)) returns the address of one-past last NumericalRadial object with angular momentum l
    const NumericalRadial** cend(const int l) const
    {
        assert(ntype_ > 0 && l >= 0 && l <= lmax_);
        return iter_ + std::accumulate(nl_, nl_ + l + 1, 0);
    }
    //!@}

    /*! @name property setters for all RadialSet objects
     *
     *  @see RadialSet
     *                                                                     */
    //!@{
    //! Set a spherical Bessel transformers for all RadialSet objects
    //! @see RadialSet::set_transformer
    void set_transformer(ModuleBase::SphericalBesselTransformer* const sbt = nullptr, const int update = 0);

    //! Set a common grid for all RadialSet objects
    //! @see RadialSet::set_grid
    void set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode = 'i');

    //! Set a common uniform grid for all RadialSet objects
    //! @see RadialSet::set_uniform_grid
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);
    //!@}

    ModuleBase::SphericalBesselTransformer* sbt() { return sbt_; }

  private:
    int ntype_ = 0;     //!< number of RadialSet in the collection
    int lmax_ = -1;     //!< maximum angular momentum of all NumericalRadial objects in the collection
    int nchi_ = 0;      //!< total number of NumericalRadial objects in the collection
    int nzeta_max_ = 0; //!< maximum number of distinct radial functions given a type & angular momentum

    //! array of RadialSet objects
    /*!
     *  Each object hold a set of NumericalRadial objects that belong to a single element.
     *  The set could be either atomic orbitals or beta functions, in which case the
     *  underlying objects are AtomicRadials or BetaRadials, respectively.
     *
     *  NOTE: AtomicRadials and BetaRadials do not necessarily have the same size as
     *  RadialSet. Therefore, a multilevel pointer is necessary for polymorphism.
     *                                                                                    */
    RadialSet** radset_ = nullptr;

    //! "Iterator" for all NumericalRadial objects
    /*!
     *   "iter_" iterates through all NumericalRadial objects from all RadialSet objects
     *   in the collection. Since NumericalRadial objects from different RadialSet objects
     *   are not contiguous, the iteration has to be done on pointers, i.e., the addresses
     *   of NumericalRadial objects are collected into a contiguous pointer array through
     *   which iter_ iterates.
     *                                                                                      */
    const NumericalRadial** iter_ = nullptr;

    //! number of NumericalRadial objects for each angular momentum
    int* nl_ = nullptr;

    //! Pointer to the object that provides spherical Bessel transforms
    /*!
     *  All NumericalRadial objects within this class should share the same
     *  spherical Bessel transformer.
     *                                                                      */
    ModuleBase::SphericalBesselTransformer* sbt_;

    //! A flag that marks the ownership of sbt_
    bool use_internal_transformer_;

    //! Deallocates all RadialSet objects and resets all members to default.
    void cleanup();

    //! Builds iter_ from radset_
    void iter_build();
};

#endif

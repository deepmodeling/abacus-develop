#ifndef RADIAL_COLLECTION_H_
#define RADIAL_COLLECTION_H_

#include <numeric>
#include <string>
#include <memory>

#include "module_basis/module_nao/radial_set.h"

/**
 * @brief A class that holds all numerical radial functions of the same kind.
 *
 * An instance of this class could be the collection of all radial functions
 * of numerical atomic orbitals, or all Kleinman-Bylander beta functions
 * from all elements involved in calculation.
 */
class RadialCollection
{
  public:
    RadialCollection() = default;
    RadialCollection(const RadialCollection& other);          ///< deep copy
    RadialCollection& operator=(const RadialCollection& rhs); ///< deep copy

    ~RadialCollection();

    /// Builds the collection from (orbital) files.
    void build(const int nfile, const std::string* const file, const char type = 'o');

    /// Builds the collection from Numerical_Nonlocal objects.
    void build(const int ntype, Numerical_Nonlocal* const nls);

    /// builds the collection from quasi hydrogen radial functions
    void build(const int ntype, const double* const charges, const int* const nmax, const std::string strategy = "minimal");

    /**
     * @name Getters
     */
    ///@{
    const std::string& symbol(const int itype) const { return radset_[itype]->symbol(); }
    int ntype() const { return ntype_; }
    int lmax(const int itype) const { return radset_[itype]->lmax(); }
    int lmax() const { return lmax_; }
    double rcut_max(const int itype) const { return radset_[itype]->rcut_max(); }
    double rcut_max() const { return rcut_max_; }
    int nzeta(const int itype, const int l) const { return radset_[itype]->nzeta(l); }
    int nzeta_max(const int itype) const { return radset_[itype]->nzeta_max(); }
    int nzeta_max() const { return nzeta_max_; }
    int nchi() const { return nchi_; }
    int nchi(const int itype) const { return radset_[itype]->nchi(); }

    const NumericalRadial& operator()(const int itype, const int l, const int izeta) const
    {
        assert(itype >= 0 && itype < ntype_);
        return radset_[itype]->chi(l, izeta);
    }

    const RadialSet& operator()(const int itype) const
    {
        assert(itype >= 0 && itype < ntype_);
        return *radset_[itype];
    }
    ///@}

    /*! @name Iterators.
     *
     *  Enable iteration through all NumericalRadial objects in the collection.
     *  Objects are sorted by l first, by itype next, by izeta last.
     */
    ///@{
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

    /// *(this->cbegin(l)) returns the address of the first NumericalRadial object with angular momentum l
    const NumericalRadial** cbegin(const int l) const
    {
        assert(ntype_ > 0 && l >= 0 && l <= lmax_);
        return iter_ + std::accumulate(nl_, nl_ + l, 0);
    }

    /// *(this->cend(l)) returns the address of one-past-last NumericalRadial object with angular momentum l
    const NumericalRadial** cend(const int l) const
    {
        assert(ntype_ > 0 && l >= 0 && l <= lmax_);
        return iter_ + std::accumulate(nl_, nl_ + l + 1, 0);
    }
    ///@}

    /**
     * @name Property setters for all RadialSet objects
     */
    ///@{
    /// Sets a spherical Bessel transformers for all RadialSet objects.
    void set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update = 0);

    /// Sets a common grid for all RadialSet objects.
    void set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode = 'i');

    /// Sets a common uniform grid for all RadialSet objects.
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);
    ///@}

  private:
    int ntype_ = 0;         ///< number of RadialSet in the collection
    int lmax_ = -1;         ///< maximum angular momentum of all NumericalRadial objects in the collection
    int nchi_ = 0;          ///< total number of NumericalRadial objects in the collection
    int nzeta_max_ = 0;     ///< maximum number of distinct radial functions given a type & angular momentum
    double rcut_max_ = 0.0; ///< maximum cutoff radius among all NumericalRadial objects

    RadialSet** radset_ = nullptr;

    /**
     * @brief "Iterator" for NumericalRadial objects.
     *
     * "iter_" iterates through all NumericalRadial objects from all RadialSet objects
     * in the collection.
     */
    const NumericalRadial** iter_ = nullptr;

    /// Number of NumericalRadial objects for each angular momentum.
    int* nl_ = nullptr;

    /// Deallocates all RadialSet objects and resets all members to default.
    void cleanup();

    /// Builds iter_ from radset_.
    void iter_build();

    /// Finds the maximum cutoff radius among all RadialSet objects and sets rcut_max_ accordingly.
    void set_rcut_max();
};

#endif

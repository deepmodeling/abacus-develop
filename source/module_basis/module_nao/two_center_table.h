#ifndef TWO_CENTER_TABLE_H
#define TWO_CENTER_TABLE_H

#include <ATen/tensor.h>
#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/radial_collection.h"

class TwoCenterTable
{
  public:
    TwoCenterTable() = default;
    ~TwoCenterTable() = default;

    TwoCenterTable(const TwoCenterTable&) = delete;
    TwoCenterTable& operator=(const TwoCenterTable&) = delete;

    void build(const RadialCollection& bra,  //!< [in] radial collection involved in <bra|op|ket>
               const RadialCollection& ket,  //!< [in] radial collection involved in <bra|op|ket>
               const char op,                //!< [in] operator of the two-center integral
               const int nr,                 //!< [in] number of table grid points
               const double cutoff           //!< [in] cutoff radius of the table
    );

    /*!
     *  @name Getters
     *                                                                                      */
    //!@{
    //! returns the operator of the two-center integral
    char op() const { return op_; }

    //! returns the number of radial points of each table
    int nr() const { return nr_; }

    // returns the number of table entries
    int ntab() const { return ntab_; }

    //! returns the radius cutoff of the table
    double rmax() const { return rmax_; }

    //! gets the read-only pointer to a specific table
    const double* table(const int itype1,        //!< [in] element index of chi1
                        const int l1,            //!< [in] angular momentum of chi1
                        const int izeta1,        //!< [in] zeta number of chi1
                        const int itype2,        //!< [in] element index of chi2
                        const int l2,            //!< [in] angular momentum of chi2
                        const int izeta2,        //!< [in] zeta number of chi2
                        const int l,             //!< [in] angular momentum of the entry
                        const bool deriv = false //!< [in] if true, return the derivative table
    ) const;

    void lookup(const int itype1,      //!< [in] element index of chi1
                const int l1,          //!< [in] angular momentum of chi1
                const int izeta1,      //!< [in] zeta number of chi1
                const int itype2,      //!< [in] element index of chi2
                const int l2,          //!< [in] angular momentum of chi2
                const int izeta2,      //!< [in] zeta number of chi2
                const int l,           //!< [in] angular momentum of the entry
                const double R,        //!< [in] distance between the two centers
                double* val,           //!< [out] interpolated values from table_
                double* dval = nullptr //!< [out] interpolated values from dtable_
    ) const;
    //!@}

    // NOTE: "nchi_ket" and "lmax_ket" are added for the purpose of "snap" in TwoCenterIntegrator
    // which generates a batch of two-center integrals depending on those information.
    // This might not be the intended purpose of this class.

    /// number of NumericalRadial objects in the ket with given itype and l
    int nchi_ket(const int itype, const int l) const
    {
        assert(itype >= 0 && itype < nchi_ket_.shape().dim_size(0));
        assert(l >= 0 && l < nchi_ket_.shape().dim_size(1));
        return nchi_ket_.get_value<int>(itype, l);
    }

    /// maximum angular momentum of the ket
    int lmax_ket() const { return nchi_ket_.shape().dim_size(1) - 1; }

  private:
    char op_ = '\0';   //!< operator associated with the present table
    int ntab_ = 0;     //!< number of table entries
    int nr_ = 0;       //!< number of radial points of each table
    double rmax_= 0.0; //!< cutoff radius of the table
    double* rgrid_ = nullptr;

    /// Table of size ntype x lmax that stores the number of radial functions of given type and l
    container::Tensor nchi_ket_{container::DataType::DT_INT, container::TensorShape({0})};

    /// two-center integral radial table, stored as a row-major matrix
    container::Tensor table_{container::DataType::DT_DOUBLE, container::TensorShape({0})};

    /// derivative table generated from cubic spline interpolation
    container::Tensor dtable_{container::DataType::DT_DOUBLE, container::TensorShape({0})};

    /// map (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index in the table
    container::Tensor index_map_{container::DataType::DT_INT, container::TensorShape({0})};

    /// returns the row-index of the table corresponding to the given two radial functions and l
    int& table_index(const NumericalRadial* ptr_rad1, const NumericalRadial* ptr_rad2, const int l);

    /// deallocates memory and reset variables to default.
    void cleanup();

    /// returns whether the given indices map to an entry in the table
    bool is_present(const int itype1,
                    const int l1,
                    const int izeta1,
                    const int itype2,
                    const int l2,
                    const int izeta2,
                    const int l) const;

    /// double factorial
    double dfact(int l) const;

    typedef void(TwoCenterTable::*looped_func)(const NumericalRadial*, const NumericalRadial*, const int l);

    /// loop-execute a function over all pairwise radial functions & l with non-vanishing Gaunt coefficients
    void two_center_loop(const RadialCollection& bra, 
                         const RadialCollection& ket, 
                         looped_func f);

    /// various looped functions during the construction of table
    void _indexing(const NumericalRadial* it1, const NumericalRadial* it2, const int l);
    void _tabulate(const NumericalRadial* it1, const NumericalRadial* it2, const int l);
};

#endif

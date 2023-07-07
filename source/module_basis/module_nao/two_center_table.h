#ifndef TWO_CENTER_TABLE_H
#define TWO_CENTER_TABLE_H

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/radial_collection.h"
#include "module_base/module_container/tensor.h"

class TwoCenterTable
{
  public:
    TwoCenterTable() {}
    ~TwoCenterTable();

    void build(const RadialCollection& bra, //!< [in] radial collection involved in <bra|op|ket>
               const RadialCollection& ket, //!< [in] radial collection involved in <bra|op|ket>
               const char op,               //!< [in] operator of the two-center integral
               const int nr,                //!< [in] number of table grid points
               const double* const rgrid,   //!< [in] table grid
               const bool deriv = false     //!< [in] whether to build the derivative of radial table
    );

    /*!
     *  @name Getters
     *                                                                                      */
    //!@{
    //! returns true if the table stores the derivative of the radial table
    bool is_deriv() const { return is_deriv_; }

    //! returns the operator of the two-center integral
    char op() const { return op_; }

    //! returns the number of radial points of each table
    int nr() const { return nr_; }

    // returns the number of table entries
    int ntab() const { return ntab_; }

    //! returns the radius cutoff of the table
    double rmax() const { return rgrid_ ? rgrid_[nr_ - 1] : 0.0; }

    //! gets the pointer to the table's grid points (read-only)
    const double* ptr_rgrid() const { return rgrid_; }

    //! gets the read-only pointer to a specific table entry
    const double* ptr_table(const int itype1, //!< [in] element index of chi1
                            const int l1,     //!< [in] angular momentum of chi1
                            const int izeta1, //!< [in] zeta number of chi1
                            const int itype2, //!< [in] element index of chi2
                            const int l2,     //!< [in] angular momentum of chi2
                            const int izeta2, //!< [in] zeta number of chi2
                            const int l       //!< [in] angular momentum of the entry
    ) const;
    //!@}

  private:
    char op_ = '\0'; //!< operator of the two-center integral
    bool is_deriv_ = false; //!< if true, table_ stores the derivative of the radial table

    int nr_ = 0; //!< number of radial points of each table
    double* rgrid_ = nullptr; //!< radial grid of each table

    int ntab_ = 0; //!< number of table entries

    //! two-center integral radial table, stored as a row-major matrix
    container::Tensor table_{ container::DataType::DT_DOUBLE, container::TensorShape({0}) };

    //! map (itype1, l1, izeta1, itype2, l2, izeta2, l) to a row index in the table
    container::Tensor index_map_{ container::DataType::DT_INT, container::TensorShape({0}) };

    //! returns the row-index of the table corresponding to the given two radial functions and l
    int& table_index(const NumericalRadial* ptr_rad1, const NumericalRadial* ptr_rad2, const int l);

    //! deallocates memory and reset variables to default.
    void cleanup();

    //! returns whether the given indices map to an entry in the table
    bool is_present(const int itype1, const int l1, const int izeta1, const int itype2, const int l2, const int izeta2, const int l) const;
};

#endif

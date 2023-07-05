#ifndef REAL_GAUNT_TABLE_H_
#define REAL_GAUNT_TABLE_H_

#include <map>
#include <array>

#include "module_base/module_container/tensor.h"

//! Table of Gaunt coefficients of real spherical harmonics
class RealGauntTable
{
public:

    RealGauntTable() {}
    ~RealGauntTable() {}

    void build(const int lmax);

    int lmax() const { return lmax_; }

    //! returns the real Gaunt coefficient
    const double& operator()(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

//private:

    int lmax_ = -1;

    //! Table of standard Gaunt coefficients
    /*!
     *  This table maps (l1,l2,l3,m1,m2,m3) to a standard Gaunt coefficient.
     *  Due to the selection rule and symmetry, only those which survive the
     *  selection rule and satisfy l1 >= l2 >= l3 && m3 >= 0 are stored.
     *                                                                                  */
    std::map<std::array<int, 6>, double> gaunt_table_;

    //! Table of real Gaunt coefficients
    /*!
     *  This table stores the real Gaunt coefficients.
     *                                                                                  */
    container::Tensor real_gaunt_table_{ container::DataType::DT_DOUBLE, container::TensorShape({0}) };

    //! Gaunt coefficients
    /*!
     *  This function computes the standard Gaunt coefficients
     *
     *                         /  m1   m2   m3
     *  G(l1,m1,l2,m2,l3,m3) = | Y    Y    Y    d Omega
     *                         /  l1   l2   l3
     *
     *  where Y is the (standard) spherical harmonics and Omega is the solid angle element.
     *
     *
     *  @note currently the algorithm evaluates the Wigner-3j symbols with the Racah
     *        formula, which might have some numerical issue for large l.
     *  @note this function computes the standard Gaunt coefficients, which is different
     *        from Gaunt coefficients of real spherical harmonics.
     *                                                                                  */
    double gaunt(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;


    //! selection rule of standard & real Gaunt coefficients about l1, l2, l3
    bool gaunt_select_l(const int l1, const int l2, const int l3) const;

    //! selection rule of standard Gaunt coefficients about m1, m2, m3
    bool gaunt_select_m(const int m1, const int m2, const int m3) const { return m1 + m2 + m3 == 0; }

    //! selection rule of real Gaunt coefficients about m1, m2, m3
    bool real_gaunt_select_m(const int m1, const int m2, const int m3) const;

    //! returns whether the given l & m are valid numbers
    /*!
     *  This function checks whether abs(mi) <= li (i=1,2,3) is satisfied.
     *  This implies li >= 0.
     *                                                                                  */
    bool is_valid_lm(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Get the Gaunt coefficient by looking up the table
    /*!
     *  This function looks up the Gaunt coefficient table. If the given parameters
     *                                                                                  */
    double gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Get the real Gaunt coefficient from the stored Gaunt coefficients
    double real_gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Symmetry-adapted key for gaunt_table_
    /*!
     *  The symmetries of Gaunt coefficients can be used to reduce the amount of storage.
     *  In particular,
     *
     *      Gaunt(l1,l2,l3,m1,m2,m3) = Gaunt(l1,l2,l3,-m1,-m2,-m3)
     *      Gaunt(1,2,3) = Gaunt(2,3,1) = Gaunt(3,1,2) = Gaunt(2,1,3) = Gaunt(1,3,2) = Gaunt(3,2,1)
     *
     *  The above symmetries enables us to store Gaunt coefficients only for
     *  l1 >= l2 >= l3 and m3 >= 0. This function permutes 1/2/3 and flips the sign
     *  of m1/m2/m3 if necessary so that the returned key (l1,l2,l3,m1,m2,m3) meets
     *  the above criteria.
     *                                                                                  */
    std::array<int, 6> gaunt_key(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! swap 1 <--> 2 if l1 < l2; do nothing otherwise
    void arrange(int& l1, int& l2, int& m1, int& m2) const;

    //! returns n! as a double
    double factorial(const int n) const;

    //! returns the linearized index of Y(l,m)
    /*!
     *  l,m    0,0   1,-1   1,0   1,1   2,-2   2,-1   2,0  ...
     *  index   0     1      2     3     4      5      6   ...
     *                                                                                  */
    int index_map(int l, int m) const;

    //! returns pow(-1, m)
    int minus_1_pow(int m) const { return m % 2 ? -1 : 1; }

};

#endif

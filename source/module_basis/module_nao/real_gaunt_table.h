#ifndef REAL_GAUNT_TABLE_H_
#define REAL_GAUNT_TABLE_H_

#include <map>
#include <array>

#include "module_base/module_container/tensor.h"

//! Table of Gaunt coefficients of real spherical harmonics
/*!
 *  This class computes and stores the Gaunt coefficients of real spherical harmonics
 *  used in two-center integrals.
 *                                                                                      */
class RealGauntTable
{
public:

    RealGauntTable() {}
    ~RealGauntTable() {}

    //! Builds the Gaunt table of real spherical harmonics
    /*!
     *  This function tabulates the Gaunt coefficients of real spherical harmonics
     *
     *                             /
     *      G(l1,l2,l3,m1,m2,m3) = | Z(l1,m1) Z(l2,m2) Z(l3,m3) d Omega
     *                             /
     *
     *  for l1,l2 <= lmax and l3 <= 2*lmax. Here Z is the real spherical harmonics
     *  defined as
     *
     *               / sqrt(2) * Re[Y(l,|m|)]   m > 0
     *               |
     *      Z(l,m) = | Y(l,0)                   m = 0
     *               |
     *               \ sqrt(2) * Im[Y(l,|m|)]   m < 0
     *
     *  @note In some literature an extra pow(-1, m) is introduced to yield a signless
     *        Cartesian expression. The definition here is consistent with
     *        ModuleBase::Ylm::sph_harm and has some minus signs, for example,
     *                       
     *             Z(1,-1) = -c * y / r
     *             Z(1, 0) = +c * z / r
     *             Z(1, 1) = -c * x / r
     *
     *        where c = sqrt(3/4/pi) and r = sqrt(x^2 + y^2 + z^2).
     *                                                                                  */
    void build(const int lmax);

    //! gets the tabulated real Gaunt coefficient
    const double& operator()(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! returns the maximum l
    int lmax() const { return lmax_; }

private:

    //! maximum angular momentum of the table (for the first two dimensions)
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
     *                             /
     *      G(l1,l2,l3,m1,m2,m3) = | Y(l1,m1) Y(l2,m2) Y(l3,m3) d Omega
     *                             /
     *
     *  where Y is the (standard) spherical harmonics and Omega is the solid angle element.
     *
     *
     *  @note This function computes the standard Gaunt coefficients, which is different
     *        from Gaunt coefficients of real spherical harmonics.
     *  @note Currently the algorithm computes the Gaunt coefficients with the Wigner-3j
     *        symbols, which in turn is evaluated with the Racah formula. This might have
     *        some numerical issue for large l and is yet to be studied later. 
     *                                                                                  */
    double gaunt(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! selection rule of standard & real Gaunt coefficients regarding l1, l2, l3
    bool gaunt_select_l(const int l1, const int l2, const int l3) const;

    //! selection rule of standard Gaunt coefficients regarding m1, m2, m3
    bool gaunt_select_m(const int m1, const int m2, const int m3) const { return m1 + m2 + m3 == 0; }

    //! selection rule of real Gaunt coefficients regarding m1, m2, m3
    bool real_gaunt_select_m(const int m1, const int m2, const int m3) const;

    //! returns whether the given l & m are valid quantum numbers
    /*!
     *  This function checks whether abs(mi) <= li (i=1,2,3) is satisfied.
     *  This implies li >= 0.
     *                                                                                  */
    bool is_valid_lm(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Get a Gaunt coefficient by looking up the table
    double gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Get a real Gaunt coefficient from the stored Gaunt coefficients
    double real_gaunt_lookup(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! Symmetry-adapted key for gaunt_table_
    /*!
     *  Standard Gaunt coefficients have the following symmetries:
     *
     *      Gaunt(l1,l2,l3,m1,m2,m3) = Gaunt(l1,l2,l3,-m1,-m2,-m3)
     *
     *      Gaunt(1,2,3) = Gaunt(2,3,1) = Gaunt(3,1,2) =
     *      Gaunt(2,1,3) = Gaunt(1,3,2) = Gaunt(3,2,1)
     *
     *  The above symmetries enable us to store merely a small portion of the Gaunt
     *  coefficients. This function permutes 1/2/3 and flips the signs of m1/m2/m3
     *  if necessary so that the returned key {l1,l2,l3,m1,m2,m3} satisfies
     *  l1 >= l2 >= l3 and m3 >= 0.
     *                                                                                  */
    std::array<int, 6> gaunt_key(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3) const;

    //! swap (l1,m1) <--> (l2,m2) if l1 < l2; do nothing otherwise
    void arrange(int& l1, int& l2, int& m1, int& m2) const;

    //! returns n! as a double
    double factorial(const int n) const;

    //! returns the linearized index of Y(l,m)
    /*!
     *  l       0   1   1   1   2   2   2   2   2   3 ...
     *  m       0  -1   0   1  -2  -1   0   1   2  -3 ...
     *  index   0   1   2   3   4   5   6   7   8   9 ...
     *                                                                                  */
    int index_map(int l, int m) const;

    //! returns pow(-1, m)
    int minus_1_pow(int m) const { return m % 2 ? -1 : 1; }

};

#endif

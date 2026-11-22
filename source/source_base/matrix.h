#ifndef MATRIX_H
#define MATRIX_H

// Peize Lin update 2018-07-29

#include <ostream>
#include<cassert>

#include <cstddef>
#include <fstream> // test
#include <limits>
#include <stdexcept>

namespace ModuleBase
{

/**
 * @brief Return nrows * ncols as size_t.
 *
 * LCAO matrices can exceed INT_MAX elements: for NLOCAL = 377884 the product
 * nrows * ncols is about 1.4e11, while INT_MAX is only 2.1e9 (issue #8342).
 * Computing the product in int overflows and turns a large allocation into a
 * negative element count, so every matrix size goes through this helper.
 */
inline std::size_t matrix_size(const int nrows, const int ncols)
{
    if (nrows < 0 || ncols < 0)
    {
        throw std::length_error("ModuleBase::matrix: negative matrix dimension");
    }

    const std::size_t rows = static_cast<std::size_t>(nrows);
    const std::size_t cols = static_cast<std::size_t>(ncols);
    if (cols != 0 && rows > std::numeric_limits<std::size_t>::max() / cols)
    {
        throw std::length_error("ModuleBase::matrix: element count overflows size_t");
    }
    return rows * cols;
}

/**
 * @brief Return the row-major index of (ir, ic) as size_t.
 *
 * The caller must guarantee that ir/ic are non-negative and in range.  The
 * multiplication is done in size_t because ir * ncols overflows int for large
 * LCAO basis sets even when the total element count still fits in int.
 */
inline std::size_t matrix_index(const int ir, const int ic, const int ncols)
{
    return static_cast<std::size_t>(ir) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(ic);
}

class matrix
{
    /* data */
public:

    int nr=0;
    int nc=0;   /* Number of rows and columns */
    double *c=nullptr;    /* Holds the data */

    /* Constructors and destructor */
    matrix(): nr(0), nc(0), c(nullptr){}
    matrix( const int nrows, const int ncols, const bool flag_zero=true );		// Peize Lin add flag_zero 2018-07-02
    matrix( const matrix &m1 ); /* copy constructor */
    matrix( matrix && m1 );			// Peize Lin add 2016-08-05
    ~matrix();

    void create( const int nrow, const int ncol, const bool flag_zero=true );			// Peize Lin add flag_zero 2018-07-02
    matrix& operator=(const matrix &m1); // Peize Lin change 2018-03-12
    matrix& operator=( matrix && m1 );	// Peize Lin add 2016-08-05

    double &operator()(const int ir,const int ic)
    {
        assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
        return c[ir*nc+ic];
    }

    const double &operator()(const int ir,const int ic) const
    {
        assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
        return c[ir*nc+ic];
    }

//	inline double &operator()(int i,int j) const
//	    { return c[nc*i+j]; }

    void operator*=(const double &s);
    void operator+=(const matrix &m);
    void operator-=(const matrix &m);

    /* member function: */
    matrix Inverse(void);

    double det(void);

    // mohan add 2011-01-13
    double trace_on(void) const;

    /* zero out all the entries */
    void zero_out(void);
    
    /* fill all entries with number */
    void fill_out(const double x);

    void get_extreme_eigen_values(double &ev_lower, double &ev_upper) const;	// mohan add 2011-01-13

    void reshape( const int nr_new, const int nc_new,  const bool flag_zero = true );		// Peize Lin add 2017-05-27

    double max() const;					// Peize Lin add 2016-09-08
    double min() const;					// Peize Lin add 2016-09-08
    double absmax() const;				// Peize Lin add 2018-07-02

    double norm() const;				// Peize Lin add 2018-08-12

    std::ostream & print( std::ostream & os, const double threshold=0.0 ) const;		// Peize Lin add 2021.09.08

    using type=double;					// Peize Lin add 2022.08.08 for template
};

matrix operator-(const matrix& m1); // unary minus
matrix operator+(const matrix &m1, const matrix &m2);
matrix operator-(const matrix &m1, const matrix &m2);
matrix operator*(const matrix &m1, const matrix &m2);
matrix operator*(const double &s, const matrix &m);
matrix operator*(const matrix &m, const double &s);

matrix transpose(const matrix &m);
double trace_on(const matrix &A, const matrix &B);		// mohan add 2011-01-13
double mdot(const matrix &A, const matrix &B);			// mohan add 2011-01-13

//std::ostream & operator<<( std::ostream & os, const matrix & m );		// Peize Lin add 2016-09-08

}

#endif // MATRIX_H

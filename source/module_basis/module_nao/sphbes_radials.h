#ifndef SPHBES_RADIALS_H_
#define SPHBES_RADIALS_H_

#include <map>
#include <vector>
#include <string>

#include "module_basis/module_nao/radial_set.h"

/**
 * @brief Numerical radials from spherical Bessel coefficients.
 *  
 */
class SphbesRadials : public RadialSet
{
  public:
    SphbesRadials() {}
    SphbesRadials(const SphbesRadials& other):
        RadialSet(other), sigma_(other.sigma_), dr_(other.dr_), coeff_(other.coeff_) {}

    SphbesRadials& operator=(const SphbesRadials& rhs);
    SphbesRadials* clone() const { return new SphbesRadials(*this); } // covariant return type

    ~SphbesRadials() {} // ~RadialSet() is called automatically

    /**
     * @brief Builds the class from a spherical Bessel coefficient file
     *
     * @param[in] file     orbital file name
     * @param[in] itype    element index in calculation
     * @param[in] dr       radial grid spacing
     * @param[in] ptr_log  output file stream for logging
     * @param[in] rank     MPI rank
     */
    void build(const std::string& file,
               const int itype = 0,
               std::ofstream* ptr_log = nullptr,
               const int rank = 0
    );

    /**
     * @brief Builds the class with truncated spherical Bessel functions.
     *
     * Instead of reading spherical Bessel coefficients from a file, an
     * alternative way to build the class is to have each truncated
     * spherical Bessel function as a NumericalRadial object.
     *
     * @param[in] lmax     maximum angular momentum
     * @param[in] nbes     number of truncated spherical Bessel functions
     * @param[in] itype    element index in calculation
     * @param[in] dr       radial grid spacing
     * @param[in] ptr_log  output file stream for logging
     * @param[in] rank     MPI rank
     */
    void build(const int lmax,
               const int nbes,
               const double rcut,
               const double sigma,
               const int itype = 0,
               std::ofstream* ptr_log = nullptr,
               const int rank = 0
    );

    void set_dr(const double dr) { dr_ = dr; }

    /**
     * @name Getters
     */
    ///@{
    double sigma() const { return sigma_; }
    double dr() const { return dr_; }
    std::vector<double> const& coeff(const int l, const int izeta) const { return coeff_.at(std::make_pair(l, izeta)); }
    ///@}

  private:

    /// Smoothing parameter.
    double sigma_ = 0.0;

    /// Radial grid spacing.
    double dr_ = 0.01;

    /// Spherical Bessel coefficients coeff_[{l,zeta}][q]
    std::map<std::pair<int,int>, std::vector<double>> coeff_;

    /// Reads spherical Bessel coefficients, cutoff radius & smoothing parameter from a file stream.
    void read_coeff(std::ifstream& ifs,
                    std::ofstream* ptr_log = nullptr,
                    const int rank = 0);

    /// 
    void build_radset(const bool normalize = true);

    /// Extracts a substring (VALUE) from a string of the form KEYWORD=" VALUE ".
    std::string extract(std::string const& str, std::string const& keyword);

    /// Splits a string into a vector of substrings with given delimiters.
    std::vector<std::string> split(std::string const& str, const char* delim = " \n\t");

    /// Computes the combination of spherical Bessel functions on a uniform grid.
    std::vector<double> sphbes_comb(const int l,
                                    std::vector<double> const& coeff_q, 
                                    double rcut,
                                    double dr,
                                    std::vector<double> const& q
    );

    /// Smoothing function.
    double smooth(double r, double rcut, double sigma);

};

#endif

